import os
from fpdf import FPDF
import requests
from collections import defaultdict
import logging
import re
import tkinter as tk
from tkinter import filedialog


# Configure logging
logging.basicConfig(level=logging.WARNING, format='%(levelname)s: %(message)s')
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')

class PathogenicDiseaseLimitReached(Exception):
    """Custom exception to indicate the limit of pathogenic diseases has been reached."""
    pass

def extract_allele(preferred_name):
    """
    Extracts the alternate allele from the preferred_name string.

    Parameters:
        preferred_name (str): The preferred name string containing allele information.

    Returns:
        str: The alternate allele if found; otherwise, "Unknown Allele".
    """
    match = re.search(r'c\.\d+([A-Z])>([A-Z])', preferred_name)
    if match:
        return match.group(2)  # Extract the alternate allele
    else:
        return "Unknown Allele"

def get_snp_info(snp_id):
    """
    Fetches ClinVar data for a given SNP ID from the myvariant.info API.

    Parameters:
        snp_id (str): The SNP ID (e.g., rs1801133).

    Returns:
        dict or list: The ClinVar data associated with the SNP ID.
    """
    url = f"https://myvariant.info/v1/variant/{snp_id}"
    params = {'fields': 'clinvar'}
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        return data
    except requests.exceptions.RequestException as e:
        logging.error(f"Error fetching data for {snp_id}: {e}")
        return None

def get_observed_alleles(data):
    """
    Extracts observed alleles from the API response.

    Parameters:
        data (dict or list): The ClinVar data fetched from the API.

    Returns:
        set: A set of observed alleles.
    """
    observed = set()
    if isinstance(data, dict):
        data = [data]
    
    for variant in data:
        clinvar_data = variant.get('clinvar', {})
        annotations = clinvar_data.get('rcv', {})
        
        if isinstance(annotations, dict):
            for entry in annotations.values():
                if isinstance(entry, dict):
                    allele = entry.get('allele') or extract_allele(entry.get('preferred_name', ''))
                    if allele and allele != "Unknown Allele":
                        observed.add(allele)
        elif isinstance(annotations, list):
            for entry in annotations:
                if isinstance(entry, dict):
                    allele = entry.get('allele') or extract_allele(entry.get('preferred_name', ''))
                    if allele and allele != "Unknown Allele":
                        observed.add(allele)
    return observed

def print_disease_info(rsid_list, observed_alleles=None, genotype=None):
    if not rsid_list or genotype is None:
        print("No disease information found.")
        return []  # Return an empty list if no data

    annotations_by_allele = defaultdict(list)
    final_pathogenic_diseases = []  # Changed to store the final pathogenic diseases found

    for variant in rsid_list:
        snp_id = variant.get('_id', 'Unknown SNP ID')  # Ensure snp_id is defined

        clinvar_data = variant.get('clinvar', {})
        annotations = clinvar_data.get('rcv', {})

        if isinstance(annotations, dict):
            annotations = [annotations]

        if isinstance(annotations, list):
            for entry in annotations:
                if not isinstance(entry, dict):
                    logging.warning(f"Skipping non-dict entry: {entry}")
                    continue

                allele = entry.get('allele') or extract_allele(entry.get('preferred_name', ''))

                # Ensure the allele matches at least one letter in the genotype
                if observed_alleles and allele not in genotype:
                    logging.debug(f"Skipping unobserved allele: {allele}")
                    continue

                accession = entry.get('accession', 'N/A')
                clinical_significance = entry.get('clinical_significance', 'Not available')

                if isinstance(clinical_significance, dict):
                    clinical_significance = clinical_significance.get('description', 'Not available')
                elif isinstance(clinical_significance, list):
                    clinical_significance = ', '.join(map(str, clinical_significance))
                elif not isinstance(clinical_significance, str):
                    clinical_significance = 'Not available'

                conditions = entry.get('conditions', [])
                condition_list = conditions if isinstance(conditions, list) else [conditions]
                disease_names = []

                for condition in condition_list:
                    disease_name = 'Unknown disease'
                    if isinstance(condition, dict):
                        name = condition.get('name') or condition.get('preferred_name')
                        if isinstance(name, str) and name.strip():
                            disease_name = name.strip()
                        else:
                            synonyms = condition.get('synonyms', [])
                            if isinstance(synonyms, list) and synonyms:
                                disease_name = ', '.join(syn.strip() for syn in synonyms if isinstance(syn, str) and syn.strip())
                    elif isinstance(condition, str):
                        disease_name = condition.strip()

                    # Ignore 'Unknown disease' and 'not provided'
                    if disease_name != 'Unknown disease' and disease_name.lower() != 'not provided':
                        disease_names.append(disease_name)

                # Check if clinical significance is 'Pathogenic' and allele matches genotype
                if clinical_significance.strip() == 'Pathogenic' and disease_names:
                    disease_names_str = ', '.join(disease_names)
                    final_pathogenic_diseases.append(disease_names_str)  # Collect the pathogenic diseases
                    
                    print(f"Found pathogenic disease for {snp_id}: {disease_names_str}")

                # Store all disease information for the allele
                annotations_by_allele[allele].append({
                    'accession': accession,
                    'disease_names': ', '.join(disease_names) if disease_names else 'not provided',
                    'clinical_significance': clinical_significance
                })

    if not annotations_by_allele:
        print("No disease annotations available.")
        return []  # Return an empty list if no annotations are found

    sorted_alleles = sorted(annotations_by_allele.keys())

    for allele in sorted_alleles:
        print(f"\nAllele: {allele}")
        print("ClinVar Accession\tDisease Names\tClinical Significance")
        for annotation in annotations_by_allele[allele]:
            print(f"{annotation['accession']}\t{annotation['disease_names']}\t{annotation['clinical_significance']}")

    # Print only the pathogenic diseases found
    print("\nList of Pathogenic Diseases:")
    if final_pathogenic_diseases:
        for disease in set(final_pathogenic_diseases):  # Use set to ensure uniqueness
            print(disease)
    else:
        print("None found.")

    return list(set(final_pathogenic_diseases))  # Return the list of unique pathogenic diseases



def process_rsids(rsid_list, genotype_list):
    all_pathogenic_diseases = []  # Cumulative list of unique pathogenic diseases

    # Filter out rsIDs that don't start with 'rs'
    rsid_list = [snp_id for snp_id in rsid_list if snp_id.startswith('rs')]

    # Create a mapping of SNP IDs to their corresponding genotypes
    snp_to_genotype = {snp_id: genotype for snp_id, genotype in zip(rsid_list, genotype_list)}

    for snp_id in rsid_list:
        genotype = snp_to_genotype.get(snp_id)
        print(f"\nProcessing SNP ID: {snp_id}, Genotype: {genotype}")

        # Fetch data for the SNP
        data = get_snp_info(snp_id)

        if data is not None:
            pathogenic_diseases = print_disease_info(data, observed_alleles=None, genotype=genotype)

            # Filter unique diseases and append them to the cumulative list
            for disease in pathogenic_diseases:
                if disease not in all_pathogenic_diseases and disease is not None:
                    all_pathogenic_diseases.append(disease)

            # Print the current pathogenic diseases found for the SNP
            if pathogenic_diseases:
                print(f"Pathogenic diseases found for {snp_id}: {', '.join(pathogenic_diseases)}")
            else:
                print(f"No pathogenic diseases found for {snp_id}.")
        else:
            print(f"Failed to retrieve data for {snp_id}.")

    # Print the final cumulative list of unique pathogenic diseases
    print("\nFinal List of Pathogenic Diseases:")
    if all_pathogenic_diseases:
        print(", ".join(all_pathogenic_diseases))
    else:
        print("None found.")

def extract_pathogenic_diseases(entry):
    # Replace this logic with how you extract diseases from each data entry
    return [disease['disease_name'] for disease in entry.get('diseases', []) if disease['clinical_significance'] == 'Pathogenic']

def open_file():
    """
    Opens a file dialog for the user to select a SNP data file and processes it.
    """
    file_path = filedialog.askopenfilename(title="Select file")
    
    if file_path:
        cumulative_pathogenic_diseases = set()  # Use a set to maintain unique diseases

        try:
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('#'):  # Skip comment lines
                        continue
                    columns = line.strip().split('\t')
                    if len(columns) > 3:  # Assuming genotype is in the fourth column
                        snp_id = columns[0]  # First column contains rsIDs
                        genotype = columns[3]  # Update this index based on your file structure

                        # Only process SNP IDs that start with "rs"
                        if not snp_id.startswith("rs"):
                            logging.info(f"Ignoring SNP ID {snp_id} as it does not start with 'rs'.")
                            continue

                        print(f"Processing SNP ID: {snp_id}, Genotype: {genotype}")
                        data = get_snp_info(snp_id)
                        if data is not None:
                            if isinstance(data, dict):
                                data = [data]

                            observed_alleles = set()
                            for entry in data:
                                observed_alleles.update(get_observed_alleles(entry))

                            pathogenic_diseases = print_disease_info(data, observed_alleles=observed_alleles, genotype=genotype)
                            
                            if pathogenic_diseases:
                                cumulative_pathogenic_diseases.update(pathogenic_diseases)  # Add found diseases to the cumulative set

                        else:
                            print(f"Failed to retrieve data for {snp_id}.")
        
        except Exception as e:
            logging.error(f"An error occurred while processing the file: {e}")

        # Print the final list of all unique pathogenic diseases found
        print("\nFinal List of Pathogenic Diseases Found Across All SNPs:")
        if cumulative_pathogenic_diseases:
            print(", ".join(cumulative_pathogenic_diseases))
        else:
            print("None found.")


def generate_pdf(pathogenic_diseases):
    """
    Generates a PDF file with specific introductory content and a list of pathogenic diseases found.

    Parameters:
        pathogenic_diseases (list): The list of pathogenic diseases to include in the PDF.
    """
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt="Results from 23&Me", ln=True, align='C')
    pdf.cell(200, 10, txt="Disclaimer: The information provided on this PDF file is not intended as medical advice.", ln=True, align='C')
    pdf.cell(200, 10, txt="Always consult with a qualified healthcare provider for medical advice and treatment.", ln=True, align='C')

    # Add a second page
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt="Top Twenty Hits", ln=True, align='C')

    # List pathogenic diseases
    pdf.set_font("Arial", size=10)
    if pathogenic_diseases:
        for disease in pathogenic_diseases:
            pdf.cell(200, 10, txt=disease, ln=True)
    else:
        pdf.cell(200, 10, txt="None found.", ln=True)

    # Save the PDF file in the Documents directory
    documents_path = os.path.expanduser("~/Documents")  # Path to the Documents folder
    pdf_output_path = os.path.join(documents_path, "ExportedResults.pdf")

    try:
        pdf.output(pdf_output_path)
        print(f"PDF containing results was exported successfully to the file named {pdf_output_path}.")
    except Exception as e:
        logging.error(f"Error saving PDF: {e}")

def open_file():
    """
    Opens a file dialog for the user to select a SNP data file and processes it.
    """
    file_path = filedialog.askopenfilename(title="Select SNP Data File")
    
    if file_path:
        cumulative_pathogenic_diseases = set()  # Use a set to maintain unique diseases
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('#'):  # Skip comment lines
                        continue
                    columns = line.strip().split('\t')
                    if len(columns) > 3:  # Assuming genotype is in the fourth column
                        snp_id = columns[0]  # First column contains rsIDs
                        genotype = columns[3]  # Update this index based on your file structure

                        # Only process SNP IDs that start with "rs"
                        if not snp_id.startswith("rs"):
                            logging.info(f"Ignoring SNP ID {snp_id} as it does not start with 'rs'.")
                            continue

                        print(f"Processing SNP ID: {snp_id}, Genotype: {genotype}")
                        data = get_snp_info(snp_id)  # Ensure this function is implemented correctly
                        if data is not None:
                            if isinstance(data, dict):
                                data = [data]

                            # Assuming this function exists and gets observed alleles
                            observed_alleles = set()
                            for entry in data:
                                observed_alleles.update(get_observed_alleles(entry))  # Ensure this function is implemented

                            pathogenic_diseases = print_disease_info(data, observed_alleles=observed_alleles, genotype=genotype)

                            if pathogenic_diseases:
                                cumulative_pathogenic_diseases.update(pathogenic_diseases)  # Add found diseases to the cumulative set

                        else:
                            print(f"Failed to retrieve data for {snp_id}.")

        except Exception as e:
            logging.error(f"An error occurred while processing the file: {e}")

        # Print the final list of all unique pathogenic diseases found
        print("\nFinal List of Pathogenic Diseases Found Across All SNPs:")
        if cumulative_pathogenic_diseases:
            print(", ".join(cumulative_pathogenic_diseases))
            generate_pdf(list(cumulative_pathogenic_diseases))  # Generate the PDF with the results
        else:
            print("None found.")

def create_ui():
    root = tk.Tk()
    root.title("SNP Raw Data Parser")
    root.configure(bg="black")

    greeting = tk.Label(root, text="SNP Raw Data Parser", font=("Arial", 20, "bold", "underline"), fg="white", bg="black")
    greeting.pack(pady=20, padx=30)

    welcome1 = tk.Label(root, text="Welcome to our SNP Raw Data Parser!", font=("Times New Roman", 11, "italic"), fg="white", bg="black")
    welcome2 = tk.Label(root, text="Please utilize the button below to upload a raw SNP data file (obtained from 23&me or ancestry.com). A report will be generated highlighting the most relevant pieces of genetic information with clinical significance.", font=("Times New Roman", 11), wraplength=500, fg="white", bg="black")
    welcome1.pack(pady=10)
    welcome2.pack()
    
    upload_button = tk.Button(root, text="Upload File", command=open_file)
    upload_button.pack(pady=20)
    
    root.mainloop()

if __name__ == "__main__":
    create_ui()


