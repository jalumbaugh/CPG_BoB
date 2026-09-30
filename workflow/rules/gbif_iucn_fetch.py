import argparse
import csv
import re
import sys
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

GBIF_SPECIES_MATCH = "https://api.gbif.org/v1/species/match"
GBIF_SPECIES = "https://api.gbif.org/v1/species/{key}"
GBIF_OCCURRENCES = "https://api.gbif.org/v1/occurrence/search"
REQUEST_TIMEOUT = (10, 120)
MAX_RETRIES = 5

CODE_TO_COUNTRY = {
    "AF": "Afghanistan", "AL": "Albania", "DZ": "Algeria", "AD": "Andorra", "AO": "Angola",
    "AR": "Argentina", "AM": "Armenia", "AU": "Australia", "AT": "Austria", "AZ": "Azerbaijan",
    "BS": "Bahamas", "BH": "Bahrain", "BD": "Bangladesh", "BB": "Barbados", "BY": "Belarus",
    "BE": "Belgium", "BZ": "Belize", "BJ": "Benin", "BT": "Bhutan", "BO": "Bolivia",
    "BA": "Bosnia and Herzegovina", "BW": "Botswana", "BR": "Brazil", "BN": "Brunei",
    "BG": "Bulgaria", "BF": "Burkina Faso", "BI": "Burundi", "KH": "Cambodia", "CM": "Cameroon",
    "CA": "Canada", "CV": "Cape Verde", "CF": "Central African Republic", "TD": "Chad",
    "CL": "Chile", "CN": "China", "CO": "Colombia", "KM": "Comoros", "CG": "Congo",
    "CR": "Costa Rica", "HR": "Croatia", "CU": "Cuba", "CY": "Cyprus", "CZ": "Czech Republic",
    "CD": "Democratic Republic of the Congo", "DK": "Denmark", "DJ": "Djibouti", "DM": "Dominica",
    "DO": "Dominican Republic", "EC": "Ecuador", "EG": "Egypt", "SV": "El Salvador",
    "GQ": "Equatorial Guinea", "ER": "Eritrea", "EE": "Estonia", "ET": "Ethiopia", "FJ": "Fiji",
    "FI": "Finland", "FR": "France", "GA": "Gabon", "GM": "Gambia", "GE": "Georgia", "DE": "Germany",
    "GH": "Ghana", "GR": "Greece", "GD": "Grenada", "GT": "Guatemala", "GN": "Guinea",
    "GW": "Guinea-Bissau", "GY": "Guyana", "HT": "Haiti", "HN": "Honduras", "HU": "Hungary",
    "IS": "Iceland", "IN": "India", "ID": "Indonesia", "IR": "Iran", "IQ": "Iraq", "IE": "Ireland",
    "IL": "Israel", "IT": "Italy", "CI": "Ivory Coast", "JM": "Jamaica", "JP": "Japan",
    "JO": "Jordan", "KZ": "Kazakhstan", "KE": "Kenya", "KI": "Kiribati", "KW": "Kuwait",
    "KG": "Kyrgyzstan", "LA": "Laos", "LV": "Latvia", "LB": "Lebanon", "LS": "Lesotho",
    "LR": "Liberia", "LY": "Libya", "LI": "Liechtenstein", "LT": "Lithuania", "LU": "Luxembourg",
    "MG": "Madagascar", "MW": "Malawi", "MY": "Malaysia", "MV": "Maldives", "ML": "Mali",
    "MT": "Malta", "MR": "Mauritania", "MU": "Mauritius", "MX": "Mexico", "MD": "Moldova",
    "MC": "Monaco", "MN": "Mongolia", "MA": "Morocco", "MZ": "Mozambique", "MM": "Myanmar",
    "NA": "Namibia", "NR": "Nauru", "NP": "Nepal", "NL": "Netherlands", "NZ": "New Zealand",
    "NI": "Nicaragua", "NE": "Niger", "NG": "Nigeria", "KP": "North Korea", "NO": "Norway",
    "OM": "Oman", "PK": "Pakistan", "PW": "Palau", "PA": "Panama", "PG": "Papua New Guinea",
    "PY": "Paraguay", "PE": "Peru", "PH": "Philippines", "PL": "Poland", "PT": "Portugal",
    "QA": "Qatar", "RO": "Romania", "RU": "Russia", "RW": "Rwanda", "KN": "Saint Kitts and Nevis",
    "LC": "Saint Lucia", "VC": "Saint Vincent and the Grenadines", "WS": "Samoa", "SM": "San Marino",
    "SA": "Saudi Arabia", "SN": "Senegal", "RS": "Serbia", "SC": "Seychelles", "SL": "Sierra Leone",
    "SG": "Singapore", "SK": "Slovakia", "SI": "Slovenia", "SB": "Solomon Islands", "SO": "Somalia",
    "ZA": "South Africa", "KR": "South Korea", "SS": "South Sudan", "ES": "Spain", "LK": "Sri Lanka",
    "SD": "Sudan", "SR": "Suriname", "SE": "Sweden", "CH": "Switzerland", "SY": "Syria",
    "TW": "Taiwan", "TJ": "Tajikistan", "TZ": "Tanzania", "TH": "Thailand", "TL": "Timor-Leste",
    "TG": "Togo", "TO": "Tonga", "TT": "Trinidad and Tobago", "TN": "Tunisia", "TR": "Turkey",
    "TM": "Turkmenistan", "TV": "Tuvalu", "UG": "Uganda", "UA": "Ukraine", "AE": "United Arab Emirates",
    "GB": "United Kingdom", "US": "United States", "UY": "Uruguay", "UZ": "Uzbekistan",
    "VU": "Vanuatu", "VE": "Venezuela", "VN": "Vietnam", "YE": "Yemen", "ZM": "Zambia", "ZW": "Zimbabwe",
    "EH": "Western Sahara", "AQ": "Antarctica"
}

session = requests.Session()
retries = Retry(
    total=MAX_RETRIES,
    backoff_factor=1,
    status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=None,
    raise_on_status=False,
)
session.mount("https://", HTTPAdapter(max_retries=retries))
session.mount("http://", HTTPAdapter(max_retries=retries))


def extract_rank_value(taxon, rank: str):
    if not isinstance(taxon, dict):
        return ""

    rank_key = rank.lower()
    for key in (rank_key, rank_key.capitalize()):
        value = taxon.get(key)
        if value:
            if isinstance(value, dict):
                return value.get("scientificName") or value.get("name") or ""
            return str(value)

    for nested_key in ("classification", "higherClassification", "higherClassificationMap"):
        nested = taxon.get(nested_key)
        if isinstance(nested, dict):
            lowered = {str(k).lower(): v for k, v in nested.items()}
            for key in (rank_key, rank_key.capitalize()):
                value = lowered.get(key)
                if value:
                    if isinstance(value, dict):
                        return value.get("scientificName") or value.get("name") or ""
                    return str(value)

    return ""


def gbif_species_lookup(scientific_name: str):
    try:
        response = session.get(
            GBIF_SPECIES_MATCH,
            params={"scientificName": scientific_name, "strict": "false"},
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()
        payload = response.json()
        key = payload.get("speciesKey") or payload.get("usageKey")
        if not key:
            return None
        taxon_response = session.get(GBIF_SPECIES.format(key=key), timeout=REQUEST_TIMEOUT)
        taxon_response.raise_for_status()
        return taxon_response.json()
    except requests.exceptions.RequestException:
        return None


def gbif_country_codes_for_taxon(taxon_key: int):
    """Fetch distinct country codes for a GBIF taxon key.

    Follows the approach described in the GBIF data-blog post on querying
    long species lists: query occurrences by taxonKey and read the country
    facet. If the facet comes back empty, fall back to paging through a
    sample of occurrence records and collecting the `country` field
    directly (still scoped to this single taxonKey, so it stays well
    within normal search-API limits without needing a bulk download).
    """
    try:
        response = session.get(
            GBIF_OCCURRENCES,
            params={
                "taxonKey": taxon_key,
                "limit": 0,
                "facet": "country",
                "facetLimit": 1000,
                "facetMinCount": 1,
            },
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()
        payload = response.json()
        for facet in payload.get("facets", []):
            if facet.get("field") in {"country", "countryCode"}:
                codes = []
                for item in facet.get("counts", []):
                    name = item.get("name")
                    if name and re.fullmatch(r"[A-Z]{2}", str(name)):
                        codes.append(str(name))
                if codes:
                    return "; ".join(sorted(set(codes)))
    except requests.exceptions.RequestException as exc:
        print(f"GBIF country facet lookup failed for taxonKey={taxon_key}: {exc}", file=sys.stderr)

    return gbif_country_codes_from_records(taxon_key)


def gbif_country_codes_from_records(taxon_key: int, max_pages: int = 3, page_size: int = 300):
    codes = set()
    try:
        for page in range(max_pages):
            response = session.get(
                GBIF_OCCURRENCES,
                params={
                    "taxonKey": taxon_key,
                    "limit": page_size,
                    "offset": page * page_size,
                },
                timeout=REQUEST_TIMEOUT,
            )
            response.raise_for_status()
            payload = response.json()
            results = payload.get("results", [])
            if not results:
                break
            for record in results:
                code = record.get("countryCode") or record.get("country")
                if code and re.fullmatch(r"[A-Z]{2}", str(code)):
                    codes.add(str(code))
            if payload.get("endOfRecords", True):
                break
    except requests.exceptions.RequestException as exc:
        print(f"GBIF record-based country lookup failed for taxonKey={taxon_key}: {exc}", file=sys.stderr)

    return "; ".join(sorted(codes))


def resolve_taxon(scientific_name: str):
    return gbif_species_lookup(scientific_name) or {}


def resolve_taxonomy(taxon):
    if not isinstance(taxon, dict):
        return "", "", ""

    class_name = extract_rank_value(taxon, "class")
    order_name = extract_rank_value(taxon, "order")
    family_name = extract_rank_value(taxon, "family")

    if not order_name and isinstance(taxon.get("props"), dict):
        class_name = extract_rank_value(taxon["props"], "class") or class_name
        order_name = extract_rank_value(taxon["props"], "order") or order_name
        family_name = extract_rank_value(taxon["props"], "family") or family_name

    return class_name, order_name, family_name


def code_list_to_names(codes):
    if not codes:
        return ""
    names = []
    for code in re.split(r"\s*;\s*", str(codes).strip()):
        if not code:
            continue
        full_name = CODE_TO_COUNTRY.get(code.upper(), code.upper())
        if full_name not in names:
            names.append(full_name)
    return "; ".join(names)


def resolve_countries(taxon: dict):
    if not isinstance(taxon, dict):
        return "", "unknown"

    key = taxon.get("key")
    if key:
        countries = gbif_country_codes_for_taxon(int(key))
        if countries:
            return code_list_to_names(countries), "gbif"

    return "", "unknown"


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def write_rows(path: Path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def parse_genus_species(row):
    cells = [cell.strip() for cell in row]
    if len(cells) < 3:
        raise SystemExit(f"Malformed input row: {row!r}")

    name_cell = cells[2]
    parts = name_cell.split(None, 1)
    if len(parts) < 2:
        raise SystemExit(f"Malformed 'genus species' value in row: {row!r}")
    genus, species = parts[0], parts[1]
    return genus, species


def enrich_taxa_file(input_path: Path, output_path: Path | None = None):
    rows = read_rows(input_path)
    if not rows:
        raise SystemExit(f"Input file is empty: {input_path}")

    # First row is always a header; the third column holds "genus species".
    data_rows = rows[1:]

    output_rows = [["Assembly Accession", "Organism Taxonomic ID", "class", "order", "family", "genus", "species", "countries", "country_source"]]
    for row in data_rows:
        if not row or all(not cell.strip() for cell in row):
            continue

        cells = [cell.strip() for cell in row]
        if len(cells) < 3:
            raise SystemExit(f"Malformed input row: {row!r}")
        assembly_accession, taxonomic_id = cells[0], cells[1]

        genus, species = parse_genus_species(row)
        scientific_name = f"{genus} {species}"
        taxon = resolve_taxon(scientific_name)

        class_name, order_name, family_name = resolve_taxonomy(taxon)
        countries, country_source = resolve_countries(taxon)

        output_rows.append([assembly_accession, taxonomic_id, class_name, order_name, family_name, genus, species, countries, country_source])

    target = output_path or input_path.with_name(f"{input_path.stem}_gbif.tsv")
    write_rows(target, output_rows)
    print(f"Wrote {len(output_rows)-1} taxa rows to {target}")


def main():
    parser = argparse.ArgumentParser(description="Add class/order/family and GBIF country names to a genus/species TSV.")
    parser.add_argument("input_tsv", help="TSV containing genus and species columns")
    parser.add_argument("--output", help="Optional output TSV path; default is <input>_gbif.tsv")
    args = parser.parse_args()

    input_path = Path(args.input_tsv).expanduser()
    if not input_path.exists():
        raise SystemExit(f"File not found: {input_path}")

    output_path = Path(args.output).expanduser() if args.output else None
    enrich_taxa_file(input_path, output_path)


# --------------------------------------------------------------------------
# Snakemake entry point
# --------------------------------------------------------------------------
# Used when this file is invoked via a Snakemake rule's `script:` directive,
# e.g.:
#
#   rule taxonomy_report:
#       input:
#           "build/.bin/acc2taxid/{project_name}_acc2taxid_lookup.tsv",
#       output:
#           "build/reports/{project_name}_taxonomy_report.tsv",
#       script:
#           "gbif_iucn_fetch.py"
#
# Snakemake execs this script with a `snakemake` object injected as a global,
# so the block below only runs under Snakemake and is skipped for normal CLI use.

def run_from_snakemake(snakemake) -> None:
    input_path = Path(snakemake.input[0]).expanduser()
    output_path = Path(snakemake.output[0]).expanduser()

    if not input_path.exists():
        raise SystemExit(f"File not found: {input_path}")

    enrich_taxa_file(input_path, output_path)


try:
    snakemake  # noqa: F821 — injected by Snakemake's `script:` directive
except NameError:
    snakemake = None

if snakemake is not None:
    run_from_snakemake(snakemake)
elif __name__ == "__main__":
    main()
