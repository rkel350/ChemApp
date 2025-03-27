import streamlit as st
from chempy import balance_stoichiometry
import re

st.title("Stoichiometry Helper App")

# ---------------------------
# Section 1: Balancing Equations
# ---------------------------
st.header("Balance Chemical Equations")

reactants_input = st.text_input("Enter reactants (comma separated):", "Fe2O3, C")
products_input = st.text_input("Enter products (comma separated):", "Fe, CO2")

if st.button("Balance Equation"):
    try:
        # Convert input strings into lists
        reactants = [r.strip() for r in reactants_input.split(',') if r.strip()]
        products = [p.strip() for p in products_input.split(',') if p.strip()]
        reac, prod = balance_stoichiometry(reactants, products)

        # Format the balanced equation output
        balanced_reactants = " + ".join(f"{reac[compound]} {compound}" for compound in reac)
        balanced_products = " + ".join(f"{prod[compound]} {compound}" for compound in prod)
        st.success(f"Balanced Equation: {balanced_reactants} → {balanced_products}")
    except Exception as e:
        st.error(f"Error balancing equation: {e}")


# ---------------------------
# Section 2: Moles-Grams Conversion
# ---------------------------

# Function to parse a chemical formula into its elemental composition.
def parse_formula(formula):
    def multiply_dict(d, factor):
        return {k: v * factor for k, v in d.items()}

    token_pattern = re.compile(r'([A-Z][a-z]?)(\d*)')
    stack = []
    current_dict = {}
    i = 0
    while i < len(formula):
        char = formula[i]
        if char == '(':
            stack.append(current_dict)
            current_dict = {}
            i += 1
        elif char == ')':
            i += 1
            num = ''
            while i < len(formula) and formula[i].isdigit():
                num += formula[i]
                i += 1
            factor = int(num) if num else 1
            current_dict = multiply_dict(current_dict, factor)
            temp = stack.pop()
            for elem, cnt in current_dict.items():
                temp[elem] = temp.get(elem, 0) + cnt
            current_dict = temp
        else:
            m = token_pattern.match(formula, i)
            if m:
                elem = m.group(1)
                cnt = int(m.group(2)) if m.group(2) else 1
                current_dict[elem] = current_dict.get(elem, 0) + cnt
                i += len(m.group(0))
            else:
                i += 1
    return current_dict


# A basic dictionary of atomic masses (in g/mol)
atomic_masses = {
    'H': 1.008, 'He': 4.0026,
    'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305, 'Al': 26.982,
    'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078,
    # Add more elements as needed
}


def compute_molar_mass(formula):
    composition = parse_formula(formula)
    total_mass = 0.0
    for element, count in composition.items():
        if element in atomic_masses:
            total_mass += atomic_masses[element] * count
        else:
            raise ValueError(f"Atomic mass for element '{element}' not found.")
    return total_mass


st.header("Moles–Grams Conversion")

formula_input = st.text_input("Enter chemical formula for conversion (e.g. H2O):", "H2O", key="formula")
conversion_type = st.selectbox("Conversion type", ["Moles to Grams", "Grams to Moles"])
amount = st.number_input("Enter amount", min_value=0.0, step=0.1, key="amount")

if st.button("Convert"):
    try:
        molar_mass = compute_molar_mass(formula_input)
        if conversion_type == "Moles to Grams":
            grams = amount * molar_mass
            st.success(
                f"{amount} moles of {formula_input} weighs {grams:.2f} grams (Molar mass: {molar_mass:.2f} g/mol).")
        else:
            moles = amount / molar_mass
            st.success(f"{amount} grams of {formula_input} is {moles:.2f} moles (Molar mass: {molar_mass:.2f} g/mol).")
    except Exception as e:
        st.error(f"Error in conversion: {e}")
