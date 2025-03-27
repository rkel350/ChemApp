import streamlit as st
from chempy import balance_stoichiometry
from chempy.util.parsing import formula_to_composition
import re

st.title("🧪 Stoichiometry Helper")

# --------------------------------------
# Section 1: Balancing Chemical Equations
# --------------------------------------

st.header("⚖️ Balance Chemical Equations")

st.markdown("""
**How to use:**
1. Type the reactants and products as chemical formulas (e.g., `Fe2O3, C` → `Fe, CO`).
2. Press **Balance Equation**.
3. Optionally check **"Show steps"** to see how it was balanced!
""")

reactants_input = st.text_input("Enter reactants (comma separated):", "Fe2O3, C")
products_input = st.text_input("Enter products (comma separated):", "Fe, CO")
show_steps = st.checkbox("Show steps")

if st.button("Balance Equation"):
    try:
        reactants = [r.strip() for r in reactants_input.split(',') if r.strip()]
        products = [p.strip() for p in products_input.split(',') if p.strip()]
        reac, prod = balance_stoichiometry(reactants, products)

        # Show balanced result
        balanced_reactants = " + ".join(f"{reac[compound]} {compound}" for compound in reac)
        balanced_products = " + ".join(f"{prod[compound]} {compound}" for compound in prod)
        st.success(f"✅ Balanced Equation: {balanced_reactants} → {balanced_products}")

        # Step-by-step breakdown
        if show_steps:
            st.markdown("### 🧠 Step-by-step Breakdown")

            def count_atoms(coeffs):
                atom_count = {}
                for compound, coeff in coeffs.items():
                    elements = formula_to_composition(compound)
                    for el, count in elements.items():
                        total = count * coeff
                        atom_count[el] = atom_count.get(el, 0) + total
                return atom_count

            reactant_atoms = count_atoms(reac)
            product_atoms = count_atoms(prod)

            st.markdown("#### 1. What the Coefficients Mean:")
            for compound, coeff in reac.items():
                st.write(f"- Reactant `{compound}` has a coefficient of {coeff}: {coeff} molecule(s).")
            for compound, coeff in prod.items():
                st.write(f"- Product `{compound}` has a coefficient of {coeff}: {coeff} molecule(s).")

            st.markdown("#### 2. Atom Counts (after applying coefficients):")
            st.write("**Reactants:**")
            for el, count in reactant_atoms.items():
                st.write(f"- {el}: {count} atom(s)")
            st.write("**Products:**")
            for el, count in product_atoms.items():
                st.write(f"- {el}: {count} atom(s)")

            st.markdown("#### 3. Balance Check:")
            balanced = True
            for el in reactant_atoms:
                r = reactant_atoms[el]
                p = product_atoms.get(el, 0)
                if r == p:
                    st.write(f"✅ {el} is balanced with {r} atoms on each side.")
                else:
                    st.write(f"❌ {el} is unbalanced: {r} in reactants vs {p} in products.")
                    balanced = False

            if balanced:
                st.success("✅ Success! All elements are balanced. Conservation of mass confirmed.")
            else:
                st.warning("⚠️ Not all elements are balanced. Double-check your formulas.")

    except Exception as e:
        st.error(f"❌ Error: {e}")

# --------------------------------------
# Section 2: Moles ↔ Grams Conversion
# --------------------------------------

st.header("⚗️ Moles ↔ Grams Conversion")

# Basic periodic table dictionary
atomic_masses = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
    'Fe': 55.845, 'Cu': 63.546, 'Zn': 65.38, 'Ag': 107.8682, 'Au': 196.966569,
    # Add more as needed
}

def parse_formula(formula):
    def multiply_dict(d, factor):
        return {k: v * factor for k, v in d.items()}

    token_pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
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

def compute_molar_mass(formula):
    composition = parse_formula(formula)
    total_mass = 0.0
    for element, count in composition.items():
        if element in atomic_masses:
            total_mass += atomic_masses[element] * count
        else:
            raise ValueError(f"Element '{element}' not found in atomic mass table.")
    return total_mass

formula_input = st.text_input("Enter a chemical formula (e.g. H2O):", "H2O")
conversion_type = st.selectbox("Conversion type:", ["Moles to Grams", "Grams to Moles"])
amount = st.number_input("Enter amount:", min_value=0.0, step=0.1)

if st.button("Convert"):
    try:
        molar_mass = compute_molar_mass(formula_input)
        if conversion_type == "Moles to Grams":
            grams = amount * molar_mass
            st.success(f"{amount} mole(s) of {formula_input} = {grams:.2f} grams (Molar mass = {molar_mass:.2f} g/mol)")
        else:
            moles = amount / molar_mass
            st.success(f"{amount} grams of {formula_input} = {moles:.4f} mole(s) (Molar mass = {molar_mass:.2f} g/mol)")
    except Exception as e:
        st.error(f"❌ Error: {e}")
