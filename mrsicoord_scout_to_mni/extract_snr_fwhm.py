import csv
import os
import re # Regular expression library

# --- Configuration ---
# 1. Path to your input CSV file
INPUT_CSV_FILE = 'points.csv' # <--- CHANGE THIS

# 2. Base path to the project directory structure
#    Set this to the part of the path BEFORE the subject/session folders
#    Example: /ix1/ginger/dsarpal/lab/reorg/projects/mrsi/Z5_mrsi/SCZ/
#    (Adjust based on your specific setup, maybe include the group folder 'SCZ' here or define it separately)
BASE_PATH = '/ix1/ginger/dsarpal/lab/reorg/projects/mrsi/Z5_mrsi/SCZ/' # <--- !!! ADJUST THIS BASE PATH !!!

# 3. Path for the output CSV file with extracted data
OUTPUT_CSV_FILE = 'extracted_mrsi_data.csv' # <--- You can change this name

# 4. Define how to construct the path to the csi.ps file
#    Based on the example: BASE_PATH/sub-{subj}/ses-{ses}/out/spectrum.{x}.{y}.dir/csi.ps
def construct_csi_ps_path(base_dir, subj, ses, x, y):
    """Constructs the full path to the target csi.ps file."""
    try:
        # Format parts of the path
        subject_folder = f"sub-{subj}"
        session_folder = f"ses-{ses}"
        spectrum_folder = f"spectrum.{x}.{y}.dir"
        csi_filename = "csi.ps"

        # Combine parts using os.path.join for cross-platform compatibility
        full_path = os.path.join(
            base_dir,
            subject_folder,
            session_folder,
            "out_old", # The 'out' directory seems constant
            spectrum_folder,
            csi_filename
        )
        return full_path
    except Exception as e:
        print(f"Error constructing path for {subj},{ses},{x},{y}: {e}")
        return None

# 5. Regular expressions based on the screenshot (LCModel output)
#    Looks for "S/N" or "FWHM", followed by optional space, '=', optional space,
#    and captures the number (integer or float).
SN_PATTERN = re.compile(r'S/N\s*=\s*([-+]?\d*\.?\d+)', re.IGNORECASE)
FWHM_PATTERN = re.compile(r'FWHM\s*=\s*([-+]?\d*\.?\d+)', re.IGNORECASE)
# --- End Configuration ---

def extract_data_from_ps(ps_filepath):
    """
    Extracts S/N and FWHM from an LCModel PostScript file using regex.

    Args:
        ps_filepath (str): The full path to the .ps file.

    Returns:
        tuple: (sn_value, fwhm_value) or (None, None) if not found or error.
    """
    sn_found = None
    fwhm_found = None

    if ps_filepath is None: # Check if path construction failed
        return None, None

    if not os.path.exists(ps_filepath):
        # This is expected sometimes if not all coordinates have a corresponding file
        # Print less verbosely or comment out if too noisy
        # print(f"Info: File not found - {ps_filepath}")
        return None, None

    try:
        # Try reading with common encodings
        encodings_to_try = ['latin-1', 'ascii', 'utf-8']
        content = None
        for enc in encodings_to_try:
            try:
                with open(ps_filepath, 'r', encoding=enc) as f:
                    content = f.read()
                break
            except UnicodeDecodeError:
                continue
            except Exception as e_inner:
                print(f"Warning: Error reading file {ps_filepath} with encoding {enc}: {e_inner}")
                continue

        if content is None:
            print(f"Error: Could not read file {ps_filepath} with tried encodings.")
            return None, None

        # Search the whole content for the patterns
        sn_match = SN_PATTERN.search(content)
        fwhm_match = FWHM_PATTERN.search(content)

        if sn_match:
            try:
                sn_found = float(sn_match.group(1))
            except (ValueError, IndexError):
                print(f"Warning: Could not parse S/N value from match '{sn_match.group(0)}' in {ps_filepath}")
                sn_found = 'parsing_error'
        # else:
            # print(f"Info: S/N pattern not found in {ps_filepath}") # Comment out if too noisy

        if fwhm_match:
            try:
                fwhm_found = float(fwhm_match.group(1))
            except (ValueError, IndexError):
                 print(f"Warning: Could not parse FWHM value from match '{fwhm_match.group(0)}' in {ps_filepath}")
                 fwhm_found = 'parsing_error'
        # else:
            # print(f"Info: FWHM pattern not found in {ps_filepath}") # Comment out if too noisy


        # Only return numeric values or None
        if isinstance(sn_found, str): sn_found = None
        if isinstance(fwhm_found, str): fwhm_found = None

        return sn_found, fwhm_found

    except FileNotFoundError:
        # Should be caught by os.path.exists, but handle defensively
        print(f"Error: File disappeared between check and open? {ps_filepath}")
        return None, None
    except Exception as e:
        print(f"Error processing file {ps_filepath}: {e}")
        return None, None

# --- Main Processing Logic ---
all_extracted_data = []
# Define header for output CSV - storing relative path/identifier now
output_header = ['SUBJ', 'SES', 'x', 'y', 'z', 't', 'label', 'spectrum_dir', 'SN', 'FWHM']

try:
    with open(INPUT_CSV_FILE, 'r', newline='') as infile:
        reader = csv.reader(infile)
        # Assuming the first row is the header
        input_header = next(reader)
        print(f"Input CSV Header: {input_header}")

        # Find column indices dynamically (more robust)
        try:
            subj_idx = input_header.index('SUBJ')
            ses_idx = input_header.index('SES')
            x_idx = input_header.index('x')
            y_idx = input_header.index('y')
            z_idx = input_header.index('z')
            t_idx = input_header.index('t')
            label_idx = input_header.index('label')
        except ValueError as e:
            print(f"Error: Missing required column in input CSV: {e}. Header found: {input_header}")
            exit()

        print(f"Processing rows from {INPUT_CSV_FILE}...")
        processed_count = 0
        found_count = 0
        for i, row in enumerate(reader):
            # Basic row validation
            if not row or len(row) <= max(subj_idx, ses_idx, x_idx, y_idx, z_idx, t_idx, label_idx):
                print(f"Warning: Skipping incomplete row {i+1}: {row}")
                continue

            try:
                # Extract identifying info from the row
                subj = row[subj_idx].strip()
                ses = row[ses_idx].strip()
                x = row[x_idx].strip()
                y = row[y_idx].strip()
                z = row[z_idx].strip()
                t = row[t_idx].strip()
                label = row[label_idx].strip()

                # Ensure coordinates are usable in path names (basic check)
                if not (subj and ses and x and y):
                    print(f"Warning: Skipping row {i+1} due to missing SUBJ/SES/X/Y: {row}")
                    continue

                # Construct the full path to the target csi.ps file
                full_ps_path = construct_csi_ps_path(BASE_PATH, subj, ses, x, y)
                processed_count += 1

                # Extract data (S/N and FWHM) from the .ps file
                sn_value, fwhm_value = extract_data_from_ps(full_ps_path)

                # Store results for this file (even if data not found, to know it was processed)
                spectrum_dir_identifier = f"spectrum.{x}.{y}.dir" if full_ps_path else "path_error"
                if sn_value is not None or fwhm_value is not None:
                    found_count += 1

                all_extracted_data.append({
                    'SUBJ': subj,
                    'SES': ses,
                    'x': x,
                    'y': y,
                    'z': z,
                    't': t,
                    'label': label,
                    'spectrum_dir': spectrum_dir_identifier, # Store the specific spectrum dir name
                    'SN': sn_value if sn_value is not None else '', # Use empty string for missing values in CSV
                    'FWHM': fwhm_value if fwhm_value is not None else ''
                })

                if processed_count % 50 == 0: # Print progress update every 50 rows
                     print(f"  ... processed {processed_count} rows")


            except Exception as e:
                 print(f"Error processing row {i+1} ({row}): {e}")
                 # Optionally append row with error status
                 all_extracted_data.append({
                    'SUBJ': row[subj_idx].strip() if len(row) > subj_idx else 'error',
                    'SES': row[ses_idx].strip() if len(row) > ses_idx else 'error',
                    'x': row[x_idx].strip() if len(row) > x_idx else 'error',
                    'y': row[y_idx].strip() if len(row) > y_idx else 'error',
                    'z': row[z_idx].strip() if len(row) > z_idx else 'error',
                    't': row[t_idx].strip() if len(row) > t_idx else 'error',
                    'label': row[label_idx].strip() if len(row) > label_idx else 'error',
                    'spectrum_dir': 'processing_error',
                    'SN': '',
                    'FWHM': ''
                 })

except FileNotFoundError:
    print(f"Error: Input CSV file not found at {INPUT_CSV_FILE}")
    exit()
except Exception as e:
    print(f"An unexpected error occurred during CSV reading or main processing loop: {e}")
    exit()

# --- Write All Results to Output CSV ---
if all_extracted_data:
    print(f"\nProcessed {processed_count} coordinates from CSV.")
    print(f"Found and extracted data for {found_count} spectrum files.")
    print(f"Writing {len(all_extracted_data)} results to {OUTPUT_CSV_FILE}...")
    try:
        with open(OUTPUT_CSV_FILE, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=output_header)
            writer.writeheader()
            writer.writerows(all_extracted_data)
        print(f"Extraction complete. Output saved to {OUTPUT_CSV_FILE}")
    except Exception as e:
        print(f"Error writing output CSV file: {e}")
else:
    print("\nNo data was processed or generated.")
