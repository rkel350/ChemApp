import streamlit as st
from chempy import balance_stoichiometry
from chempy.util.parsing import formula_to_composition
import re
from PIL import Image
import random



if "balanced" not in st.session_state:
    st.session_state.balanced = False
logo = Image.open("bitschemlabpurp.png")
st.image(logo, width=200)

st.title("🧪 Stoichiometry Helper")

# --------------------------------------
# Section 1: Balancing Chemical Equations
# --------------------------------------

st.header("⚖️ Balance Chemical Equations")

st.markdown("""
**Enter an unbalanced equation using this format:**  
👉 `Fe2O3 + C -> Fe + CO`
""")

equation_input = st.text_input("Enter unbalanced equation:", "Fe2O3 + C -> Fe + CO")
show_steps = st.checkbox("Show steps")
show_bca = st.checkbox("Show B-C-A Table (Before - Change - After)")


# Try to split the equation into reactants and products
try:
    if "->" in equation_input:
        left_side, right_side = equation_input.split("->")
        reactants = [r.strip() for r in left_side.split("+")]
        products = [p.strip() for p in right_side.split("+")]
    else:
        reactants = products = []
        st.warning("Please use '->' to separate reactants and products.")
except Exception as e:
    st.error(f"❌ Could not parse equation: {e}")
    reactants = products = []

# Only show the balance button if parsing succeeded
if reactants and products:
    if st.button("Balance Equation"):
        st.session_state.balanced = True
    if st.session_state.balanced:
        try:
            reac, prod = balance_stoichiometry(reactants, products)

            # Show balanced result
            balanced_reactants = " + ".join(f"{reac[compound]} {compound}" for compound in reac)
            balanced_products = " + ".join(f"{prod[compound]} {compound}" for compound in prod)
            st.success(f"✅ Balanced Equation: {balanced_reactants} → {balanced_products}")

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

                st.markdown("#### 1. Coefficients Explained:")
                for compound, coeff in reac.items():
                    st.write(f"- Reactant `{compound}`: {coeff} molecule(s)")
                for compound, coeff in prod.items():
                    st.write(f"- Product `{compound}`: {coeff} molecule(s)")

                st.markdown("#### 2. Atom Counts After Balancing:")
                st.write("**Reactants:**")
                for el, count in reactant_atoms.items():
                    st.write(f"- {el}: {count} atom(s)")
                st.write("**Products:**")
                for el, count in product_atoms.items():
                    st.write(f"- {el}: {count} atom(s)")

                st.markdown("#### 3. Final Check:")
                balanced = True
                for el in reactant_atoms:
                    r = reactant_atoms[el]
                    p = product_atoms.get(el, 0)
                    if r == p:
                        st.write(f"✅ {el} is balanced with {r} atoms on both sides.")
                    else:
                        st.write(f"❌ {el} is unbalanced: {r} in reactants vs {p} in products.")
                        balanced = False

                if balanced:
                    st.success("✅ All atoms are balanced!")
                else:
                    st.warning("⚠️ Not all atoms are balanced. Double-check the equation.")


            if show_bca:
                st.markdown("### 🧪 B-C-A Table Setup")

                # Let user input starting moles for each substance
                starting_moles = {}
                st.markdown("#### Enter starting moles for each substance:")
                all_substances = set(list(reac.keys()) + list(prod.keys()))
                for compound in all_substances:
                    starting_moles[compound] = st.number_input(
                        f"{compound} (mol)", min_value=0.0, value=0.0, step=0.1
                    )

                # Find limiting reactant based on stoichiometric ratios
                limiting_ratios = {
                    compound: starting_moles[compound] / coeff
                    for compound, coeff in reac.items()
                    if coeff > 0 and starting_moles[compound] > 0
                }

                if not limiting_ratios:
                    st.warning("Please enter starting moles for at least one reactant.")
                else:
                    limiting_compound = min(limiting_ratios, key=limiting_ratios.get)
                    reaction_scale = limiting_ratios[limiting_compound]
                    st.info(
                        f"🔍 Limiting reactant: **{limiting_compound}** (based on {reaction_scale:.2f} reaction scale)")

                    # Build BCA Table
                    st.markdown("### 📊 B-C-A Table")

                    import pandas as pd

                    data = {"Substance": [], "Before (mol)": [], "Change (mol)": [], "After (mol)": []}

                    for compound in all_substances:
                        coeff_reac = reac.get(compound, 0)
                        coeff_prod = prod.get(compound, 0)
                        coeff = coeff_prod - coeff_reac  # Net change

                        before = starting_moles.get(compound, 0.0)
                        change = coeff * reaction_scale
                        after = before + change

                        data["Substance"].append(compound)
                        data["Before (mol)"].append(before)
                        data["Change (mol)"].append(change)
                        data["After (mol)"].append(after)

                    df_bca = pd.DataFrame(data)
                    st.dataframe(df_bca.style.format(precision=2))

        except Exception as e:
            st.error(f"❌ Error balancing equation: {e}")
else:
    st.info("Waiting for a properly formatted equation...")

# --------------------------------------
# Section 2: Moles ↔ Grams Conversion
# --------------------------------------

st.header("⚗️ Moles ↔ Grams Conversion")

atomic_masses = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
    'Fe': 55.845, 'Cu': 63.546, 'Zn': 65.38, 'Ag': 107.8682, 'Au': 196.966569,
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
            raise ValueError(f"Element '{element}' not found.")
    return total_mass

formula_input = st.text_input("Enter a chemical formula (e.g. H2O):", "H2O")
conversion_type = st.selectbox("Conversion type:", ["Moles to Grams", "Grams to Moles"])
amount = st.number_input("Enter amount:", min_value=0.0, step=0.1)

if st.button("Convert"):
    try:
        molar_mass = compute_molar_mass(formula_input)
        if conversion_type == "Moles to Grams":
            grams = amount * molar_mass
            st.success(f"{amount} mole(s) of {formula_input} = {grams:.2f} grams")
        else:
            moles = amount / molar_mass
            st.success(f"{amount} grams of {formula_input} = {moles:.4f} mole(s)")
    except Exception as e:
        st.error(f"❌ Error: {e}")

# Random derpy footer lines
footers = [
    "👀 If you got this far, you probably deserve a snack.",
    "🚫 No goggles were harmed in the making of this app.",
    "🔬 Built with 90% caffeine and 10% panic.",
    "🔥 If your eyebrows are still intact, congrats!",
    "🧼 Wash your hands. You touched sulfur.",
    "💥 Hit reload to simulate an explosion.",
    "🧪 Slightly more stable than the average group project.",
    "⚗️ Available in mole-sized servings.",
    "🥽 No safety goggles? No problem... just kidding. Wear them."
]

random_footer = random.choice(footers)

# Show the footer
st.markdown(f"""
<hr style="margin-top: 3em; margin-bottom: 1em;">
<div style="text-align: center; font-size: 16px; color: gray;">
    <p>{random_footer}</p>
    <p style="font-size: 12px;">© Bit’s Chem Lab, 2025. All unstable ions reserved.</p>
</div>
""", unsafe_allow_html=True)


