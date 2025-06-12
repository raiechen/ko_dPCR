# Import the necessary libraries
import streamlit as st
import pandas as pd
import io # Use io to handle the uploaded file bytes
import numpy as np
from io import BytesIO

def load_data(uploaded_file):
    """
    Reads an uploaded CSV file, handles BOM, custom separators, and returns a DataFrame.
    It tries to read with 'utf-8-sig' and falls back to 'latin-1' on error.
    """
    file_buffer = BytesIO(uploaded_file.getvalue())
    
    try:
        # Try to read with utf-8-sig first
        first_line_str = file_buffer.readline().decode('utf-8-sig').strip()
        file_buffer.seek(0)  # Reset after peeking

        if first_line_str.startswith("sep="):
            delimiter = first_line_str.split("=")[1]
            
            df = pd.read_csv(file_buffer, sep=delimiter, skiprows=1, header=0, encoding='utf-8-sig')
        else:
            
            df = pd.read_csv(file_buffer, header=0, encoding='utf-8-sig')
        
        return df

    except UnicodeDecodeError:
        
        file_buffer.seek(0)  # Reset buffer before trying again

        # Fallback to latin-1
        first_line_str = file_buffer.readline().decode('latin-1').strip()
        file_buffer.seek(0)
        
        if first_line_str.startswith("sep="):
            delimiter = first_line_str.split("=")[1]
            
            df = pd.read_csv(file_buffer, sep=delimiter, skiprows=1, header=0, encoding='latin-1')
        else:
            
            df = pd.read_csv(file_buffer, header=0, encoding='latin-1')
        
        return df

#start with v71 in box: Benchling_KO_dPCR_auto_analysis.py
#criteia variables
Max_CI_selected_target = 25
Max_CI_B2M =10
LLOQ_selected_target = 3
NTC_pos_part_selected_target = 10
NTC_pos_part_B2M = 10
Min_valid_part_selected_target = 40000
Min_valid_part_B2M = 40000
FS_non_edit_max_cutoff = 10
Ref_lot_non_edit_max = 105
Ref_lot_non_edit_min = 95 



# Set the title of the Streamlit app
st.title("TRAC/TRBC2/TGFBR2/TRBC2-9 KO dPCR Analysis App")

# Add a dropdown for target selection
target_options = ['Select a Target','TGFBR2-LNA', 'TRAC-LNA', 'TRBC2-LNA', 'TRBC2-9-LNA']
selected_target = st.selectbox("Select the Target for Analysis:", target_options)

# Check if the user has selected a valid target
if selected_target == 'Select a Target':
    st.warning("Please select a valid target for analysis.")
else:
    # Proceed with the analysis only if a valid target is selected
    st.write(f"You selected: {selected_target}")
    # Add your analysis logic here

# Add a file uploader widget
uploaded_file = st.file_uploader("Choose a CSV file", type=['csv'])



# Check if a file has been uploaded
if uploaded_file is not None:
    try:
        df = load_data(uploaded_file)

        
        if df is None or df.empty:
            st.warning("The uploaded file was loaded, but the DataFrame is empty.")
        else:
            # Always display the uploaded DataFrame in a collapsible section
            with st.expander("Show Uploaded DataFrame"):
                st.dataframe(df)

        # --- Handling for new CSV format ---
        # Define the column mapping from the new format to the expected format
        column_mapping = {
            'Target (Name)': 'Target',
            'Conc. [cp/µL] (dPCR reaction)': 'Conc. [copies/µL]',
            'CI (95%) (dPCR reaction)': 'CI (95%)',
            'Partitions (Valid)': 'Partitions (valid)',
            'Partitions (Positive)': 'Partitions (positive)'
        }
        # Check if any of the new format's columns are present
        if any(col in df.columns for col in column_mapping.keys()):
            
            df.rename(columns=column_mapping, inplace=True)
            # In the new format, some data is spread across rows. Forward fill the 'Sample/NTC/Control' column.
            df['Sample/NTC/Control'] = df['Sample/NTC/Control'].ffill()

        # Ensure required columns exist
        required_columns = ['Sample/NTC/Control', 'Target', 'Conc. [copies/µL]', 'CI (95%)']
        if not all(col in df.columns for col in required_columns):
            st.error(f"Missing required columns in the uploaded file: {required_columns}")
            print(f"Missing required columns in the uploaded file: {required_columns}")
        else:
            # Filter data for the selected target and 'B2M'
            target_data = df[df['Target'] == selected_target]
            b2m_data = df[df['Target'] == 'B2M']

            # Merge the data on 'Sample/NTC/Control' to align the selected target and B2M rows
            merged_data = pd.merge(
                target_data[['Sample/NTC/Control', 'Conc. [copies/µL]','CI (95%)','Partitions (valid)']],
                b2m_data[['Sample/NTC/Control', 'Conc. [copies/µL]', 'CI (95%)','Partitions (valid)']],
                on='Sample/NTC/Control',
                suffixes=(f'_{selected_target}', '_B2M')
            )

            # Calculate %WT
            merged_data['%WT'] = (
                (merged_data[f'Conc. [copies/µL]_{selected_target}'] /
                 merged_data['Conc. [copies/µL]_B2M']) * 100
            ).round(2)

            
            # Create new_df with relevant columns
            new_df = merged_data[['Sample/NTC/Control', '%WT',f'Conc. [copies/µL]_{selected_target}',f'Conc. [copies/µL]_B2M', f'CI (95%)_{selected_target}', f'CI (95%)_B2M',f'Partitions (valid)_{selected_target}', 'Partitions (valid)_B2M']].rename(
                columns={
                    'Sample/NTC/Control': 'Sample Description',
                    f'Conc. [copies/µL]_{selected_target}': 'Target Concentration (copies/µL)',
                    f'CI (95%)_{selected_target}': 'Target CI (95%)',
                    f'Conc. [copies/µL]_B2M': 'Ref Concentration (copies/µL)',
                    f'CI (95%)_B2M': 'Ref CI (95%)',
                    f'Partitions (valid)_{selected_target}': 'Target Partitions (valid)',
                    'Partitions (valid)_B2M': 'Ref Partitions (valid)'
                }
            )
            new_df['Target CI (95%)'] = new_df['Target CI (95%)'].str.replace('%', '', regex=False)
            new_df['Ref CI (95%)'] = new_df['Ref CI (95%)'].str.replace('%', '', regex=False)

            # Add a new column '%KO' = 100 - '%WT'
            new_df['%KO'] = (100 - new_df['%WT']).round(2)
            # Filter out rows where 'Sample Description' is 'NTC'
            new_df = new_df[new_df['Sample Description'] != 'NTC']
            new_df = new_df[new_df['Sample Description'].str.lower() != 'pdna only']

            # Define conditions for 'Designation'
            # Define conditions for 'Designation'
            conditions = [
                # Condition for 'D'
                (
                    (new_df['Target CI (95%)'].astype(float) <= Max_CI_selected_target) &
                    (new_df['Ref CI (95%)'].astype(float) <= Max_CI_B2M) &
                    (new_df['Target Concentration (copies/µL)'] >= LLOQ_selected_target) &
                    (new_df['Target Partitions (valid)'] >= Min_valid_part_selected_target) &
                    (new_df['Ref Partitions (valid)'] >= Min_valid_part_B2M) &
                    (~np.isinf(new_df['%WT']))  # Ensure %WT is not infinity

                    
                ),
                # Condition for 'NQ'
                (
                    (new_df['Ref CI (95%)'].astype(float) <= Max_CI_B2M) &
                    (new_df['Target Concentration (copies/µL)'] < LLOQ_selected_target) &
                    (new_df['Target Partitions (valid)'] >= Min_valid_part_selected_target) &
                    (new_df['Ref Partitions (valid)'] >= Min_valid_part_B2M) &
                    (~np.isinf(new_df['%WT']))  # Ensure %WT is not infinity

                ),
                # Condition for 'UND'
                (
                    (new_df['Ref CI (95%)'].astype(float) > Max_CI_B2M) |
                    (new_df['Target Partitions (valid)'] < Min_valid_part_selected_target) |
                    (new_df['Ref Partitions (valid)'] < Min_valid_part_B2M) |
                    (np.isinf(new_df['%WT']))  # Add condition for %WT being infinity

                )
            ]

            # Define corresponding values for each condition
            values = ['D', 'NQ', 'UND']
            # Add the 'Designation' column to new_df
            new_df['Designation'] = np.select(conditions, values, default='')

            # Reorder columns to place '%KO' right after '%WT'
            columns = list(new_df.columns)  # Get the current column order
            wt_index = columns.index('%WT')  # Find the index of the '%WT' column
            # Insert '%KO' right after '%WT'
            columns.insert(wt_index + 1, columns.pop(columns.index('%KO')))
            new_df = new_df[columns]  # Reorder the DataFrame
            


            # Display a success message
            st.success(f"File '{uploaded_file.name}' successfully processed and analyzed!")
            # Check assay acceptance criteria for NTC
            ntc_data = df[df['Sample/NTC/Control'] == 'NTC']  # Filter rows where 'Sample/NTC/Control' is 'NTC'

            if not ntc_data.empty:
                # Check if 'Partitions (positive)' <= 10 for all NTC rows
                if ntc_data['Partitions (positive)'].max() <= 10:
                    # Assay passes
                    st.info("Assay Status: Pass", icon="✅")  # Green highlight
                else:
                    # Assay fails
                    st.error("Assay Status: Fail", icon="❌")  # Red highlight
            else:
                st.warning("No NTC data found in the uploaded file.")
            
            st.write("Criteria used for analysis:")
            st.write(f"Max CI for selected target: {Max_CI_selected_target}")
            st.write(f"Max CI for B2M: {Max_CI_B2M}")
            st.write(f"LLOQ for selected target: {LLOQ_selected_target}")
            st.write(f"NTC positive partitions for selected target: {NTC_pos_part_selected_target}")
            st.write(f"NTC positive partitions for B2M: {NTC_pos_part_B2M}")
            st.write(f"Min valid partitions for selected target: {Min_valid_part_selected_target}")
            st.write(f"Min valid partitions for B2M: {Min_valid_part_B2M}")
            st.write(f"FS non-edit max cutoff: {FS_non_edit_max_cutoff}")
            st.write(f"Ref lot non-edit max: {Ref_lot_non_edit_max}")
            st.write(f"Ref lot non-edit min: {Ref_lot_non_edit_min}")
            st.write("---")
            st.write("Analyzed Data:")
            st.dataframe(new_df)
            # Add a button to export the data as an Excel report
            if st.button("Export as Excel Report"):
                try:
                    # Add assay status to the report
                    ntc_data = df[df['Sample/NTC/Control'] == 'NTC']  # Filter rows where 'Sample/NTC/Control' is 'NTC'
                    if not ntc_data.empty and ntc_data['Partitions (positive)'].max() <= 10:
                        assay_status = "Pass"
                    else:
                        assay_status = "Fail"

                    # Add assay status as a new row in the report
                    report_df = new_df.copy()
                    report_df.loc[len(report_df)] = ["Assay Status", assay_status] + [""] * (len(report_df.columns) - 2)

                    # Save the DataFrame to an Excel file in memory
                    output = BytesIO()
                    with pd.ExcelWriter(output, engine='openpyxl') as writer:  # Use openpyxl as the engine
                        report_df.to_excel(writer, index=False, sheet_name="Analysis Report")
                        

                    # Provide the Excel file for download
                    st.download_button(
                        label="Download Excel Report",
                        data=output.getvalue(),
                        file_name="Analysis_Report.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                except Exception as e:
                    st.error(f"Error exporting Excel report: {e}")
                    print(f"Error exporting Excel report: {e}")

        
    

    except Exception as e:
        # Display a more general error message if something else goes wrong
        st.error(f"Error processing file: {e}")
        # Also print error to terminal for debugging
        print(f"Error processing file: {e}")

else:
    # Display a message prompting the user to upload a file in the web app
    st.info("Please upload a CSV file.")