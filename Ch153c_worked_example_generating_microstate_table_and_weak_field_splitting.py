#%% CONFIGURATION SECTION - MODIFY THESE VALUES
# L defines the orbital angular momentum (0=s, 1=p, 2=d, 3=f, etc.)
L = 2  
# r defines the number of electrons
r = 3 
################################################################################

#%% IMPORTS AND SETUP
import csv
import os
from itertools import product, combinations as combo_iter
from collections import defaultdict, Counter

# Orbital labels corresponding to L values
ORBITAL_LABELS = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N", "O", "Q", "R", "T", "U", "V", "W", "X", "Y", "Z"]

# Define the possible values for m_l and m_s based on L
ML_VALUES = list(range(-L, L+1))       # Possible m_l values for each electron
SPIN_VALUES = [+0.5, -0.5]             # Possible spin states (m_s)
NUM_ELECTRONS = r                      # Number of electrons

print(f"Configuration: L={L} ({ORBITAL_LABELS[L]} orbital), {r} electrons")
print(f"Possible m_l values: {ML_VALUES}")
print(f"Possible m_s values: {SPIN_VALUES}")

#%% GENERATING MICROSTATES
# Generate all possible combinations of m_l and m_s for each electron
all_electron_states = list(product(ML_VALUES, SPIN_VALUES))

# Get unique combinations of electron states (order doesn't matter)
unique_state_combinations = list(combo_iter(all_electron_states, r))

# Convert to the flat format used in the rest of the code
filtered_combinations = []
for state_combo in unique_state_combinations:
    flat_combo = [val for electron in state_combo for val in electron]
    filtered_combinations.append(flat_combo)

# Format a state using the requested convention
def format_state(combo):
    state = []
    for i in range(r):
        ml = combo[2*i]
        ms = combo[2*i+1]
        spin_symbol = "+" if ms == 0.5 else "-"
        state.append(f"{ml}{spin_symbol}")
    return "(" + ", ".join(state) + ")"

# Group states by M_L and M_S
grouped_states = defaultdict(list)

for combo in filtered_combinations:
    M_L = sum(combo[2*i] for i in range(r))
    M_S = sum(combo[2*i+1] for i in range(r))
    formatted_state = format_state(combo)
    grouped_states[(M_L, M_S)].append(formatted_state)

#%% GETTING TERM SYMBOLS
def get_term_symbols(l, r):
    """Return a list of term symbols for the configuration l^r."""

    # Total number of (ml, ms) pairs for this subshell.
    n = (2*l+1)*2
    # All possible values of ml = -l, -l+1, ..., l-1, l.
    ml = list(range(-l,l+1))
    # All possible values of 2ms = -1, 1. That is, ms = -1/2, +1/2. We work
    # with 2ms instead of ms so that we can handle integers only.
    ms2 = [-1,1]
    # All possible (ml, 2ms) pairs for this subshell.
    ml_ms2 = list(product(ml, ms2))

    # All possible microstates for r electrons in this subshell.
    microstates = list(combo_iter(range(n), r))
    # The totals ML = sum(ml) and MS2 = sum(2ms) for each microstate
    ML = [sum([ml_ms2[microstate[j]][0] for j in range(r)])
                                    for microstate in microstates]
    MS2 = [sum([ml_ms2[microstate[j]][1] for j in range(r)])
                                    for microstate in microstates]
    # Count the microstates (MS, ML). Store them this way round so we can
    # pick off the ground state term (maximum S) first.
    MS2_ML = Counter(zip(MS2,ML))
    N = len(microstates)

    # Extract the term symbols by starting at the minimum (ML, MS) value and
    # removing microstates corresponding to the (L, S) term it belongs to.
    # Repeat until we're out of microstates.
    terms = []
    while N>0:
        S, L = min(MS2_ML)
        terms.append('{}{}'.format(-S+1, ORBITAL_LABELS[-L]))
        for ML in range(L, -L+1):
            for MS in range(S, -S+1,2):
                MS2_ML[MS,ML] -= 1
                if MS2_ML[MS,ML] == 0:
                    del MS2_ML[MS,ML]
                N -= 1
    return terms

# Get term symbols for the current configuration
term_symbols_list = get_term_symbols(L, r)
print(f"Term symbols for {ORBITAL_LABELS[L]}^{r}: {', '.join(term_symbols_list)}")

#%% FUNCTION TO CHECK IF A CELL BELONGS TO A TERM SYMBOL
def is_term_symbol(ml, ms, multiplicity, orbital_letter):
    """
    Check if a cell with given M_L and M_S belongs to a specific term symbol.
    
    Args:
        ml: The M_L value of the cell
        ms: The M_S value of the cell
        multiplicity: The multiplicity (2S+1) of the term
        orbital_letter: The orbital letter (S, P, D, F, G, H, etc.)
    
    Returns:
        True if the cell belongs to the specified term symbol, False otherwise
    """
    # Calculate S from multiplicity
    S = (multiplicity - 1) / 2
    
    # Calculate L from orbital letter
    L = ORBITAL_LABELS.index(orbital_letter)
    
    # Check if M_S is in the valid range for this S value
    valid_ms_values = [i/2 for i in range(-int(2*S), int(2*S)+1, 1)]
    if ms not in valid_ms_values:
        return False
    
    # Check if M_L is in the valid range for this L value
    valid_ml_values = list(range(-L, L+1))
    if ml not in valid_ml_values:
        return False
    
    return True

#%% DICTIONARY OF TERM SYMBOLS WITH COLORS
term_symbols = {
    # S orbital terms
    "1S": {"multiplicity": 1, "orbital": "S", "color": "#FF0000"},  # Red
    "2S": {"multiplicity": 2, "orbital": "S", "color": "#0000FF"},  # Blue
    "3S": {"multiplicity": 3, "orbital": "S", "color": "#00FF00"},  # Green
    "4S": {"multiplicity": 4, "orbital": "S", "color": "#FF00FF"},  # Magenta
    "5S": {"multiplicity": 5, "orbital": "S", "color": "#800000"},  # Maroon
    "6S": {"multiplicity": 6, "orbital": "S", "color": "#008000"},  # Dark Green
    "7S": {"multiplicity": 7, "orbital": "S", "color": "#800080"},  # Purple
    "8S": {"multiplicity": 8, "orbital": "S", "color": "#FFA500"},  # Orange
    
    # P orbital terms
    "1P": {"multiplicity": 1, "orbital": "P", "color": "#800080"},  # Purple
    "2P": {"multiplicity": 2, "orbital": "P", "color": "#FFA500"},  # Orange
    "3P": {"multiplicity": 3, "orbital": "P", "color": "#00FFFF"},  # Cyan
    "4P": {"multiplicity": 4, "orbital": "P", "color": "#008080"},  # Teal
    "5P": {"multiplicity": 5, "orbital": "P", "color": "#FF69B4"},  # Hot Pink
    "6P": {"multiplicity": 6, "orbital": "P", "color": "#4B0082"},  # Indigo
    "7P": {"multiplicity": 7, "orbital": "P", "color": "#A52A2A"},  # Brown
    "8P": {"multiplicity": 8, "orbital": "P", "color": "#808000"},  # Olive
    
    # D orbital terms
    "1D": {"multiplicity": 1, "orbital": "D", "color": "#FF00FF"},  # Magenta
    "2D": {"multiplicity": 2, "orbital": "D", "color": "#A52A2A"},  # Brown
    "3D": {"multiplicity": 3, "orbital": "D", "color": "#808000"},  # Olive
    "4D": {"multiplicity": 4, "orbital": "D", "color": "#008080"},  # Teal
    "5D": {"multiplicity": 5, "orbital": "D", "color": "#000080"},  # Navy
    "6D": {"multiplicity": 6, "orbital": "D", "color": "#FF6347"},  # Tomato
    "7D": {"multiplicity": 7, "orbital": "D", "color": "#32CD32"},  # Lime Green
    "8D": {"multiplicity": 8, "orbital": "D", "color": "#800000"},  # Maroon
    
    # F orbital terms
    "1F": {"multiplicity": 1, "orbital": "F", "color": "#800000"},  # Maroon
    "2F": {"multiplicity": 2, "orbital": "F", "color": "#FFD700"},  # Gold
    "3F": {"multiplicity": 3, "orbital": "F", "color": "#9370DB"},  # Medium Purple
    "4F": {"multiplicity": 4, "orbital": "F", "color": "#3CB371"},  # Medium Sea Green
    "5F": {"multiplicity": 5, "orbital": "F", "color": "#CD853F"},  # Peru
    "6F": {"multiplicity": 6, "orbital": "F", "color": "#4682B4"},  # Steel Blue
    "7F": {"multiplicity": 7, "orbital": "F", "color": "#D2691E"},  # Chocolate
    "8F": {"multiplicity": 8, "orbital": "F", "color": "#8A2BE2"},  # Blue Violet
    
    # G orbital terms
    "1G": {"multiplicity": 1, "orbital": "G", "color": "#2E8B57"},  # Sea Green
    "2G": {"multiplicity": 2, "orbital": "G", "color": "#DAA520"},  # Goldenrod
    "3G": {"multiplicity": 3, "orbital": "G", "color": "#4169E1"},  # Royal Blue
    "4G": {"multiplicity": 4, "orbital": "G", "color": "#8B4513"},  # Saddle Brown
    "5G": {"multiplicity": 5, "orbital": "G", "color": "#FF4500"},  # Orange Red
    "6G": {"multiplicity": 6, "orbital": "G", "color": "#2F4F4F"},  # Dark Slate Gray
    "7G": {"multiplicity": 7, "orbital": "G", "color": "#FF8C00"},  # Dark Orange
    "8G": {"multiplicity": 8, "orbital": "G", "color": "#BDB76B"},  # Dark Khaki
    
    # H orbital terms
    "1H": {"multiplicity": 1, "orbital": "H", "color": "#556B2F"},  # Dark Olive Green
    "2H": {"multiplicity": 2, "orbital": "H", "color": "#B8860B"},  # Dark Goldenrod
    "3H": {"multiplicity": 3, "orbital": "H", "color": "#483D8B"},  # Dark Slate Blue
    "4H": {"multiplicity": 4, "orbital": "H", "color": "#8B0000"},  # Dark Red
    "5H": {"multiplicity": 5, "orbital": "H", "color": "#008B8B"},  # Dark Cyan
    "6H": {"multiplicity": 6, "orbital": "H", "color": "#9932CC"},  # Dark Orchid
    "7H": {"multiplicity": 7, "orbital": "H", "color": "#8B008B"},  # Dark Magenta
    "8H": {"multiplicity": 8, "orbital": "H", "color": "#E9967A"},  # Dark Salmon
    
    # I orbital terms
    "1I": {"multiplicity": 1, "orbital": "I", "color": "#00CED1"},  # Dark Turquoise
    "2I": {"multiplicity": 2, "orbital": "I", "color": "#FF1493"},  # Deep Pink
    "3I": {"multiplicity": 3, "orbital": "I", "color": "#00FA9A"},  # Medium Spring Green
    "4I": {"multiplicity": 4, "orbital": "I", "color": "#1E90FF"},  # Dodger Blue
    "5I": {"multiplicity": 5, "orbital": "I", "color": "#B22222"},  # Fire Brick
    "6I": {"multiplicity": 6, "orbital": "I", "color": "#228B22"},  # Forest Green
    "7I": {"multiplicity": 7, "orbital": "I", "color": "#4B0082"},  # Indigo
    "8I": {"multiplicity": 8, "orbital": "I", "color": "#20B2AA"},  # Light Sea Green
    
    # K orbital terms
    "1K": {"multiplicity": 1, "orbital": "K", "color": "#ADFF2F"},  # Green Yellow
    "2K": {"multiplicity": 2, "orbital": "K", "color": "#CD5C5C"},  # Indian Red
    "3K": {"multiplicity": 3, "orbital": "K", "color": "#F08080"},  # Light Coral
    "4K": {"multiplicity": 4, "orbital": "K", "color": "#90EE90"},  # Light Green
    "5K": {"multiplicity": 5, "orbital": "K", "color": "#FFB6C1"},  # Light Pink
    "6K": {"multiplicity": 6, "orbital": "K", "color": "#FFA07A"},  # Light Salmon
    "7K": {"multiplicity": 7, "orbital": "K", "color": "#DB7093"},  # Pale Violet Red
    "8K": {"multiplicity": 8, "orbital": "K", "color": "#FFDAB9"},  # Peach Puff
    
    # L orbital terms
    "1L": {"multiplicity": 1, "orbital": "L", "color": "#4169E1"},  # Royal Blue
    "2L": {"multiplicity": 2, "orbital": "L", "color": "#8B4513"},  # Saddle Brown
    "3L": {"multiplicity": 3, "orbital": "L", "color": "#FA8072"},  # Salmon
    "4L": {"multiplicity": 4, "orbital": "L", "color": "#F4A460"},  # Sandy Brown
    "5L": {"multiplicity": 5, "orbital": "L", "color": "#2E8B57"},  # Sea Green
    "6L": {"multiplicity": 6, "orbital": "L", "color": "#A0522D"},  # Sienna
    "7L": {"multiplicity": 7, "orbital": "L", "color": "#87CEEB"},  # Sky Blue
    "8L": {"multiplicity": 8, "orbital": "L", "color": "#6A5ACD"},  # Slate Blue
    
    # M orbital terms
    "1M": {"multiplicity": 1, "orbital": "M", "color": "#D8BFD8"},  # Thistle
    "2M": {"multiplicity": 2, "orbital": "M", "color": "#FF6347"},  # Tomato
    "3M": {"multiplicity": 3, "orbital": "M", "color": "#40E0D0"},  # Turquoise
    "4M": {"multiplicity": 4, "orbital": "M", "color": "#EE82EE"},  # Violet
    "5M": {"multiplicity": 5, "orbital": "M", "color": "#F5DEB3"},  # Wheat
    "6M": {"multiplicity": 6, "orbital": "M", "color": "#FFFF00"},  # Yellow
    "7M": {"multiplicity": 7, "orbital": "M", "color": "#9ACD32"},  # Yellow Green
    "8M": {"multiplicity": 8, "orbital": "M", "color": "#6B8E23"},  # Olive Drab
    
    # N orbital terms
    "1N": {"multiplicity": 1, "orbital": "N", "color": "#DDA0DD"},  # Plum
    "2N": {"multiplicity": 2, "orbital": "N", "color": "#B0E0E6"},  # Powder Blue
    "3N": {"multiplicity": 3, "orbital": "N", "color": "#800080"},  # Purple
    "4N": {"multiplicity": 4, "orbital": "N", "color": "#663399"},  # Rebecca Purple
    "5N": {"multiplicity": 5, "orbital": "N", "color": "#FF0000"},  # Red
    "6N": {"multiplicity": 6, "orbital": "N", "color": "#BC8F8F"},  # Rosy Brown
    "7N": {"multiplicity": 7, "orbital": "N", "color": "#4169E1"},  # Royal Blue
    "8N": {"multiplicity": 8, "orbital": "N", "color": "#8B4513"},  # Saddle Brown
    
    # O orbital terms
    "1O": {"multiplicity": 1, "orbital": "O", "color": "#6A5ACD"},  # Slate Blue
    "2O": {"multiplicity": 2, "orbital": "O", "color": "#708090"},  # Slate Gray
    "3O": {"multiplicity": 3, "orbital": "O", "color": "#FFFAFA"},  # Snow
    "4O": {"multiplicity": 4, "orbital": "O", "color": "#00FF7F"},  # Spring Green
    "5O": {"multiplicity": 5, "orbital": "O", "color": "#4682B4"},  # Steel Blue
    "6O": {"multiplicity": 6, "orbital": "O", "color": "#D2B48C"},  # Tan
    "7O": {"multiplicity": 7, "orbital": "O", "color": "#008080"},  # Teal
    "8O": {"multiplicity": 8, "orbital": "O", "color": "#D8BFD8"},  # Thistle
    
    # Q orbital terms
    "1Q": {"multiplicity": 1, "orbital": "Q", "color": "#FF6347"},  # Tomato
    "2Q": {"multiplicity": 2, "orbital": "Q", "color": "#40E0D0"},  # Turquoise
    "3Q": {"multiplicity": 3, "orbital": "Q", "color": "#EE82EE"},  # Violet
    "4Q": {"multiplicity": 4, "orbital": "Q", "color": "#F5DEB3"},  # Wheat
    "5Q": {"multiplicity": 5, "orbital": "Q", "color": "#FFFF00"},  # Yellow
    "6Q": {"multiplicity": 6, "orbital": "Q", "color": "#9ACD32"},  # Yellow Green
    "7Q": {"multiplicity": 7, "orbital": "Q", "color": "#FF0000"},  # Red
    "8Q": {"multiplicity": 8, "orbital": "Q", "color": "#0000FF"},  # Blue
    
    # R orbital terms
    "1R": {"multiplicity": 1, "orbital": "R", "color": "#00FF00"},  # Green
    "2R": {"multiplicity": 2, "orbital": "R", "color": "#FF00FF"},  # Magenta
    "3R": {"multiplicity": 3, "orbital": "R", "color": "#800000"},  # Maroon
    "4R": {"multiplicity": 4, "orbital": "R", "color": "#008000"},  # Dark Green
    "5R": {"multiplicity": 5, "orbital": "R", "color": "#800080"},  # Purple
    "6R": {"multiplicity": 6, "orbital": "R", "color": "#FFA500"},  # Orange
    "7R": {"multiplicity": 7, "orbital": "R", "color": "#00FFFF"},  # Cyan
    "8R": {"multiplicity": 8, "orbital": "R", "color": "#008080"},  # Teal

    # T orbital terms
    "1T": {"multiplicity": 1, "orbital": "T", "color": "#FF69B4"},  # Hot Pink
    "2T": {"multiplicity": 2, "orbital": "T", "color": "#4B0082"},  # Indigo
    "3T": {"multiplicity": 3, "orbital": "T", "color": "#FF00FF"},  # Magenta
    "4T": {"multiplicity": 4, "orbital": "T", "color": "#A52A2A"},  # Brown
    "5T": {"multiplicity": 5, "orbital": "T", "color": "#808000"},  # Olive
    "6T": {"multiplicity": 6, "orbital": "T", "color": "#008080"},  # Teal
    "7T": {"multiplicity": 7, "orbital": "T", "color": "#000080"},  # Navy
    "8T": {"multiplicity": 8, "orbital": "T", "color": "#FF6347"},  # Tomato
    
    # U orbital terms
    "1U": {"multiplicity": 1, "orbital": "U", "color": "#32CD32"},  # Lime Green
    "2U": {"multiplicity": 2, "orbital": "U", "color": "#800000"},  # Maroon
    "3U": {"multiplicity": 3, "orbital": "U", "color": "#FFD700"},  # Gold
    "4U": {"multiplicity": 4, "orbital": "U", "color": "#9370DB"},  # Medium Purple
    "5U": {"multiplicity": 5, "orbital": "U", "color": "#3CB371"},  # Medium Sea Green
    "6U": {"multiplicity": 6, "orbital": "U", "color": "#CD853F"},  # Peru
    "7U": {"multiplicity": 7, "orbital": "U", "color": "#4682B4"},  # Steel Blue
    "8U": {"multiplicity": 8, "orbital": "U", "color": "#D2691E"},  # Chocolate
    
    # V orbital terms
    "1V": {"multiplicity": 1, "orbital": "V", "color": "#8A2BE2"},  # Blue Violet
    "2V": {"multiplicity": 2, "orbital": "V", "color": "#2E8B57"},  # Sea Green
    "3V": {"multiplicity": 3, "orbital": "V", "color": "#DAA520"},  # Goldenrod
    "4V": {"multiplicity": 4, "orbital": "V", "color": "#4169E1"},  # Royal Blue
    "5V": {"multiplicity": 5, "orbital": "V", "color": "#8B4513"},  # Saddle Brown
    "6V": {"multiplicity": 6, "orbital": "V", "color": "#FF4500"},  # Orange Red
    "7V": {"multiplicity": 7, "orbital": "V", "color": "#2F4F4F"},  # Dark Slate Gray
    "8V": {"multiplicity": 8, "orbital": "V", "color": "#FF8C00"},  # Dark Orange
    
    # W orbital terms
    "1W": {"multiplicity": 1, "orbital": "W", "color": "#BDB76B"},  # Dark Khaki
    "2W": {"multiplicity": 2, "orbital": "W", "color": "#6495ED"},  # Cornflower Blue
    "3W": {"multiplicity": 3, "orbital": "W", "color": "#556B2F"},  # Dark Olive Green
    "4W": {"multiplicity": 4, "orbital": "W", "color": "#B8860B"},  # Dark Goldenrod
    "5W": {"multiplicity": 5, "orbital": "W", "color": "#483D8B"},  # Dark Slate Blue
    "6W": {"multiplicity": 6, "orbital": "W", "color": "#8B0000"},  # Dark Red
    "7W": {"multiplicity": 7, "orbital": "W", "color": "#008B8B"},  # Dark Cyan
    "8W": {"multiplicity": 8, "orbital": "W", "color": "#9932CC"},  # Dark Orchid
    
    # X orbital terms
    "1X": {"multiplicity": 1, "orbital": "X", "color": "#8B008B"},  # Dark Magenta
    "2X": {"multiplicity": 2, "orbital": "X", "color": "#E9967A"},  # Dark Salmon
    "3X": {"multiplicity": 3, "orbital": "X", "color": "#F0E68C"},  # Khaki
    "4X": {"multiplicity": 4, "orbital": "X", "color": "#BC8F8F"},  # Rosy Brown
    "5X": {"multiplicity": 5, "orbital": "X", "color": "#00CED1"},  # Dark Turquoise
    "6X": {"multiplicity": 6, "orbital": "X", "color": "#FF1493"},  # Deep Pink
    "7X": {"multiplicity": 7, "orbital": "X", "color": "#00FA9A"},  # Medium Spring Green
    "8X": {"multiplicity": 8, "orbital": "X", "color": "#1E90FF"},  # Dodger Blue
    
    # Y orbital terms
    "1Y": {"multiplicity": 1, "orbital": "Y", "color": "#B22222"},  # Fire Brick
    "2Y": {"multiplicity": 2, "orbital": "Y", "color": "#228B22"},  # Forest Green
    "3Y": {"multiplicity": 3, "orbital": "Y", "color": "#4B0082"},  # Indigo
    "4Y": {"multiplicity": 4, "orbital": "Y", "color": "#20B2AA"},  # Light Sea Green
    "5Y": {"multiplicity": 5, "orbital": "Y", "color": "#87CEEB"},  # Sky Blue
    "6Y": {"multiplicity": 6, "orbital": "Y", "color": "#778899"},  # Light Slate Gray
    "7Y": {"multiplicity": 7, "orbital": "Y", "color": "#B0C4DE"},  # Light Steel Blue
    "8Y": {"multiplicity": 8, "orbital": "Y", "color": "#ADFF2F"},  # Green Yellow
}

#%% CREATING OUTPUT DIRECTORY
# Create output directory with configuration name
output_dir = f"{ORBITAL_LABELS[L]}{r}_microstates"
os.makedirs(output_dir, exist_ok=True)
print(f"Output directory: {os.path.abspath(output_dir)}")

#%% WRITING BASIC MICROSTATES TABLE
# Output file for the table
output_filename = os.path.join(output_dir, 'microstates_table.csv')

# Write to CSV
with open(output_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write header
    writer.writerow(['M_L', 'M_S', 'Number of States', 'States'])
    
    # Write each group
    for (M_L, M_S), states in sorted(grouped_states.items()):
        writer.writerow([M_L, M_S, len(states), "\n".join(states)])

print(f"Microstates table saved to '{output_filename}'.")
print(f"Total possible electron states: {len(all_electron_states)}")
print(f"Number of unique microstates: {len(filtered_combinations)}")
print(f"Number of unique (M_L, M_S) groups: {len(grouped_states)}")

#%% CREATING VISUAL TABLE
# Find the range of M_L and M_S values
ml_values = sorted(set(ml for ml, _ in grouped_states.keys()))
ms_values = sorted(set(ms for _, ms in grouped_states.keys()), reverse=True)

# Create a visual table representation - TRANSPOSED
visual_table_filename = os.path.join(output_dir, 'microstates_visual_table.csv')

with open(visual_table_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write header row with M_S values
    header = ['M_L / M_S'] + [str(ms) for ms in ms_values]
    writer.writerow(header)
    
    # Write each row for different M_L values
    for ml in ml_values:
        row = [str(ml)]
        for ms in ms_values:
            if (ml, ms) in grouped_states:
                states = grouped_states[(ml, ms)]
                cell_content = f"{len(states)} states"
                row.append(cell_content)
            else:
                row.append('')
        writer.writerow(row)

print(f"Visual microstates table saved to '{visual_table_filename}'.")

#%% CREATING DETAILED TABLE
# Create a detailed table with all states listed - TRANSPOSED
detailed_table_filename = os.path.join(output_dir, 'microstates_detailed_table.csv')

with open(detailed_table_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write header row with M_S values
    header = ['M_L / M_S'] + [str(ms) for ms in ms_values]
    writer.writerow(header)
    
    # Write each row for different M_L values
    for ml in ml_values:
        row = [str(ml)]
        for ms in ms_values:
            if (ml, ms) in grouped_states:
                states = grouped_states[(ml, ms)]
                cell_content = "\n".join(states)
                row.append(cell_content)
            else:
                row.append('')
        writer.writerow(row)

print(f"Detailed microstates table saved to '{detailed_table_filename}'.")

#%% CREATING DETAILED TABLE
# Create a detailed table with all states listed - TRANSPOSED
detailed_table_filename = os.path.join(output_dir, 'microstates_detailed_table.csv')

with open(detailed_table_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write header row with M_S values
    header = ['M_L / M_S'] + [str(ms) for ms in ms_values]
    writer.writerow(header)
    
    # Write each row for different M_L values
    for ml in ml_values:
        row = [str(ml)]
        for ms in ms_values:
            if (ml, ms) in grouped_states:
                states = grouped_states[(ml, ms)]
                cell_content = "\n".join(states)
                row.append(cell_content)
            else:
                row.append('')
        writer.writerow(row)

print(f"Detailed microstates table saved to '{detailed_table_filename}'.")

# Create an HTML version of the detailed table
detailed_html_filename = os.path.join(output_dir, 'microstates_detailed_table.html')

with open(detailed_html_filename, 'w') as file:
    # Start HTML document
    file.write(f'''<!DOCTYPE html>
<html>
<head>
    <title>Detailed Microstates Table for {ORBITAL_LABELS[L]}^{NUM_ELECTRONS}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            font-size: 14px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin-bottom: 30px;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: center;
            vertical-align: top;
        }}
        th {{
            background-color: #f2f2f2;
            font-weight: bold;
        }}
        .state {{
            margin: 2px 0;
            font-family: monospace;
        }}
        .config-info {{
            margin-bottom: 15px;
            font-weight: bold;
        }}
        .ml-column {{
            font-weight: bold;
            background-color: #f8f8f8;
        }}
    </style>
</head>
<body>
    <h2>Detailed Microstates Table for {ORBITAL_LABELS[L]}^{NUM_ELECTRONS}</h2>
    
    <div class="config-info">
        Configuration: {ORBITAL_LABELS[L]}^{NUM_ELECTRONS} (L={L}, {NUM_ELECTRONS} electrons)
    </div>
    
    <table>
        <tr>
            <th>M_L / M_S</th>''')
    
    # Write header row with M_S values
    for ms in ms_values:
        file.write(f'\n            <th>{ms}</th>')
    file.write('\n        </tr>')
    
    # Write each row for different M_L values
    for ml in ml_values:
        file.write(f'\n        <tr>\n            <td class="ml-column">{ml}</td>')
        
        for ms in ms_values:
            if (ml, ms) in grouped_states:
                states = grouped_states[(ml, ms)]
                file.write('\n            <td>')
                for state in states:
                    file.write(f'\n                <div class="state">{state}</div>')
                file.write('\n            </td>')
            else:
                file.write('\n            <td></td>')
        
        file.write('\n        </tr>')
    
    # End HTML document
    file.write('''
    </table>
    
    <p>Each cell contains all possible microstates with the given M_L and M_S values.</p>
    <p>Format: (ml+/-, ml+/-, ...) where ml is the magnetic quantum number and +/- indicates spin up/down.</p>
</body>
</html>''')

print(f"Detailed microstates HTML table saved to '{detailed_html_filename}'.")



#%% CREATING HTML TABLES FOR SPECIFIC TERM SYMBOLS
def create_html_table(filename, term_symbol_to_highlight):
    """
    Create an HTML table with highlighting for a specific term symbol.
    
    Args:
        filename: The output filename
        term_symbol_to_highlight: The term symbol to highlight (e.g., "1P")
    """
    term_info = term_symbols.get(term_symbol_to_highlight)
    if not term_info:
        print(f"Term symbol {term_symbol_to_highlight} not found in the dictionary.")
        return
    
    # Get full path for the output file
    full_path = os.path.abspath(filename)
    
    with open(filename, 'w') as file:
        # Start HTML document
        file.write(f'''<!DOCTYPE html>
<html>
<head>
    <title>Microstates Table with Highlighted {term_symbol_to_highlight} Term Symbol</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
        }}
        table {{
            border-collapse: separate;
            border-spacing: 4px;
            width: 100%;
            margin-bottom: 30px;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: center;
            position: relative;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        .highlighted {{
            outline: 3px solid {term_info["color"]};
            outline-offset: 2px;
            z-index: 1;
        }}
        .legend {{
            display: flex;
            flex-wrap: wrap;
            margin-bottom: 20px;
            padding: 10px;
            border: 1px solid #ddd;
            background-color: #f9f9f9;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            margin-right: 20px;
        }}
        .legend-color {{
            width: 20px;
            height: 20px;
            margin-right: 5px;
            border: 2px solid {term_info["color"]};
        }}
        .config-info {{
            margin-bottom: 15px;
            font-weight: bold;
        }}
    </style>
</head>
<body>
    <h2>Microstates Table with Highlighted {term_symbol_to_highlight} Term Symbol</h2>
    
    <div class="config-info">
        Configuration: {ORBITAL_LABELS[L]}^{r} ({L=}, {r} electrons)
    </div>
    
    <div class="legend">
        <h3>Legend:</h3>
        <div class="legend-item">
            <div class="legend-color"></div>
            <div>{term_symbol_to_highlight} ({term_info["multiplicity"]}{term_info["orbital"]})</div>
        </div>
    </div>
    
    <table>''')
        
        # Write header row
        file.write('\n        <tr>\n')
        file.write(f'            <th>M_L / M_S</th>\n')
        for ms in ms_values:
            file.write(f'            <th>{ms}</th>\n')
        file.write('        </tr>\n')
        
        # Write data rows
        for ml in ml_values:
            file.write('        <tr>\n')
            file.write(f'            <th>{ml}</th>\n')
            
            for ms in ms_values:
                # Check if this cell belongs to the specified term symbol
                is_highlighted = is_term_symbol(ml, ms, term_info["multiplicity"], term_info["orbital"])
                highlight_class = ' class="highlighted"' if is_highlighted else ''
                
                if (ml, ms) in grouped_states:
                    states = grouped_states[(ml, ms)]
                    cell_content = f"{len(states)} states"
                    file.write(f'            <td{highlight_class}>{cell_content}</td>\n')
                else:
                    file.write(f'            <td{highlight_class}></td>\n')
            
            file.write('        </tr>\n')
        
        # End HTML document
        file.write('''    </table>
    
    <h3>Explanation:</h3>
    <p>Each cell in the table represents a specific (M_L, M_S) combination. The colored border around a cell indicates that it belongs to the term symbol shown in the legend.</p>
    <p>The number inside each cell shows how many microstates exist with that specific (M_L, M_S) combination.</p>
</body>
</html>''')

    print(f"HTML table with highlighted {term_symbol_to_highlight} term symbol saved to: {full_path}")
    return full_path

#%% CREATING HTML TABLE WITH ALL TERM SYMBOLS
def create_all_terms_html_table(filename):
    """
    Create an HTML table with all term symbols highlighted with different colored borders,
    along with a legend explaining which color corresponds to which term symbol.
    """
    # Get full path for the output file
    full_path = os.path.abspath(filename)
    
    # Filter term_symbols to only include those that are relevant for this configuration
    relevant_term_symbols = {term: info for term, info in term_symbols.items() 
                            if term in term_symbols_list}
    
    with open(filename, 'w') as file:
        # Start HTML document
        file.write('''<!DOCTYPE html>
<html>
<head>
    <title>Microstates Table with All Term Symbols</title>
    <style>
        body {
            font-family: Arial, sans-serif;
        }
        table.data-table {
            border-collapse: separate;
            border-spacing: 4px;
            width: 100%;
            margin-bottom: 30px;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: center;
            position: relative;
        }
        th {
            background-color: #f2f2f2;
        }
        .legend {
            display: flex;
            flex-wrap: wrap;
            margin-bottom: 20px;
            padding: 10px;
            border: 1px solid #ddd;
            background-color: #f9f9f9;
        }
        .legend-item {
            display: flex;
            align-items: center;
            margin-right: 20px;
            margin-bottom: 10px;
        }
        .legend-color {
            width: 20px;
            height: 20px;
            margin-right: 5px;
            border: 2px solid;
        }
        .config-info {
            margin-bottom: 15px;
            font-weight: bold;
        }
''')
        
        # Add styles for each term symbol's border
        border_styles = []
        for i, (term, info) in enumerate(relevant_term_symbols.items()):
            # Calculate offset for stacking borders
            offset = i * 2
            border_styles.append(f'''
        .{term.lower()}-border {{
            outline: 2px solid {info["color"]};
            outline-offset: {offset}px;
            z-index: {i+1};
        }}''')
        
        file.write(''.join(border_styles))
        
        file.write('''
    </style>
</head>
<body>
    <h2>Microstates Table with All Term Symbols</h2>
    
    <div class="config-info">
        Configuration: ''')
        file.write(f"{ORBITAL_LABELS[L]}^{r} (L={L}, {r} electrons)")
        file.write('''
    </div>
    
    <div class="legend">
        <h3>Legend:</h3>
        <div style="display: flex; flex-wrap: wrap; width: 100%;">''')
        
        # Add legend items only for relevant term symbols
        for term, info in relevant_term_symbols.items():
            file.write(f'''
            <div class="legend-item">
                <div class="legend-color" style="border-color: {info["color"]};"></div>
                <div>{term} ({info["multiplicity"]}{info["orbital"]})</div>
            </div>''')
        
        file.write('''
        </div>
    </div>
    
    <table class="data-table">''')
        
        # Write header row
        file.write('\n        <tr>\n')
        file.write('            <th>M_L / M_S</th>\n')
        for ms in ms_values:
            file.write(f'            <th>{ms}</th>\n')
        file.write('        </tr>\n')
        
        # Write data rows
        for ml in ml_values:
            file.write('        <tr>\n')
            file.write(f'            <th>{ml}</th>\n')
            
            for ms in ms_values:
                # Determine which term symbols this cell belongs to
                cell_classes = []
                
                for term, info in relevant_term_symbols.items():
                    if is_term_symbol(ml, ms, info["multiplicity"], info["orbital"]):
                        cell_classes.append(f"{term.lower()}-border")
                
                class_attr = f' class="{" ".join(cell_classes)}"' if cell_classes else ''
                
                if (ml, ms) in grouped_states:
                    states = grouped_states[(ml, ms)]
                    cell_content = f"{len(states)} states"
                    file.write(f'            <td{class_attr}>{cell_content}</td>\n')
                else:
                    file.write(f'            <td{class_attr}></td>\n')
            
            file.write('        </tr>\n')
        
        # End HTML document
        file.write('''    </table>
    
    <h3>Explanation:</h3>
    <p>Each cell in the table represents a specific (M_L, M_S) combination. The colored borders around a cell indicate which term symbols that cell belongs to, according to the legend above.</p>
    <p>The number inside each cell shows how many microstates exist with that specific (M_L, M_S) combination.</p>
</body>
</html>''')

    print(f"HTML table with all term symbols (overlaying borders) saved to: {full_path}")
    return full_path

#%% GENERATE HTML FILES
# Create HTML tables for specific term symbols
html_files = []

# Only create HTML files for term symbols that are relevant to this configuration
relevant_term_symbols = [term for term in term_symbols_list if term in term_symbols]

for term in relevant_term_symbols:
    filename = os.path.join(output_dir, f"microstates_{term}_term.html")
    html_path = create_html_table(filename, term)
    html_files.append(html_path)

# Create HTML table with all term symbols
all_terms_filename = os.path.join(output_dir, "microstates_all_terms.html")
all_terms_path = create_all_terms_html_table(all_terms_filename)
html_files.append(all_terms_path)

#%% PRINT SUMMARY
# Print summary of all files created
print("\nSummary of all files created:")
print(f"1. Microstates table: {os.path.abspath(output_filename)}")
print(f"2. Visual microstates table: {os.path.abspath(visual_table_filename)}")
print(f"3. Detailed microstates table: {os.path.abspath(detailed_table_filename)}")
print("4. HTML files:")
for i, file_path in enumerate(html_files):
    print(f"   {i+1}. {file_path}")

print(f"\nAll files are saved in: {os.path.abspath(output_dir)}")

###########################################################################################################

#%% Change this section only
import numpy as np

import numpy as np

def character_identity(L_value):
    """χ(E) = 2L + 1"""
    return 2 * L_value + 1

def character_Cn(L_value, phi):
    """χ[C(φ)] = sin((L + 1/2)φ) / sin(φ / 2)"""
    return np.sin((L_value + 0.5) * phi) / np.sin(phi / 2)

def character_inversion(L_value):
    """χ(i) = ±(2L + 1)"""
    return (-1)**(L*r) * (2 * L_value + 1)

def character_Sn(L_value, phi):
    """χ[S(φ)] = ± sin((L + 1/2)(φ + π)) / sin((φ + π)/2)"""
    return (-1)**(L*r) * np.sin((L_value + 0.5) * (phi + np.pi)) / np.sin((phi + np.pi) / 2)

def character_sigma(L_value):
    """χ(σ) = P + sin((L + 1/2)π) = P (-1)^L"""
    return (-1)**(L*r) * np.sin((L_value + 0.5) * np.pi)
############################################################################################################

def calculate_oh_reducible_rep(L_value):
    """Calculate the reducible representation for the Oh point group"""
    return {
        "E": character_identity(L_value),
        "8C3": character_Cn(L_value, 2*np.pi/3),
        "6C2": character_Cn(L_value, np.pi),
        "6C4": character_Cn(L_value, np.pi/2),
        "3C2'": character_Cn(L_value, np.pi),
        "i": character_inversion(L_value),
        "8S6": character_Sn(L_value, np.pi/3),
        "6σd": character_sigma(L_value),
        "6S4": character_Sn(L_value, np.pi/2),
        "3σh": character_sigma(L_value)
    }

def calculate_d4h_reducible_rep(L_value):
    """Calculate the reducible representation for the D4h point group"""
    return {
        "E": character_identity(L_value),
        "2C4": character_Cn(L_value, np.pi/2),
        "C2": character_Cn(L_value, np.pi),
        "2C2'": character_Cn(L_value, np.pi),
        "2C2''": character_Cn(L_value, np.pi),
        "i": character_inversion(L_value),
        "2S4": character_Sn(L_value, np.pi/2),
        "σh": character_sigma(L_value),
        "2σv": character_sigma(L_value),
        "2σd": character_sigma(L_value)
    }

def calculate_d5d_reducible_rep(L_value):
    """Calculate the reducible representation for the D5d point group"""
    return {
        "E": character_identity(L_value),
        "2C5": character_Cn(L_value, 2*np.pi/5),
        "2C5^2": character_Cn(L_value, 2*2*np.pi/5),
        "5C2": character_Cn(L_value, np.pi),
        "i": character_inversion(L_value),
        "2S10": character_Sn(L_value, np.pi/5),
        "2S10^3": character_Sn(L_value, 3*2*np.pi/10),
        "5σd": character_sigma(L_value)
    }


def decompose_into_irreps(reducible_rep, point_group):
    """Decompose a reducible representation into irreducible representations"""
    if point_group == "Oh":
        # Character table for Oh point group
        character_table = {
            "A1g": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            "A2g": [1, 1, -1, -1, 1, 1, 1, -1, -1, 1],
            "Eg": [2, -1, 0, 0, 2, 2, -1, 0, 0, 2],
            "T1g": [3, 0, -1, 1, -1, 3, 0, -1, 1, -1],
            "T2g": [3, 0, 1, -1, -1, 3, 0, 1, -1, -1],
            "A1u": [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
            "A2u": [1, 1, -1, -1, 1, -1, -1, 1, 1, -1],
            "Eu": [2, -1, 0, 0, 2, -2, 1, 0, 0, -2],
            "T1u": [3, 0, -1, 1, -1, -3, 0, 1, -1, 1],
            "T2u": [3, 0, 1, -1, -1, -3, 0, -1, 1, 1]
        }
        class_sizes = [1, 8, 6, 6, 3, 1, 8, 6, 6, 3]
        group_order = 48
        classes = ["E", "8C3", "6C2", "6C4", "3C2'", "i", "8S6", "6σd", "6S4", "3σh"]
    
    elif point_group == "D4h":
        # Character table for D4h point group
        character_table = {
            "A1g": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            "A2g": [1, 1, 1, -1, -1, 1, 1, 1, -1, -1],
            "B1g": [1, -1, 1, 1, -1, 1, -1, 1, 1, -1],
            "B2g": [1, -1, 1, -1, 1, 1, -1, 1, -1, 1],
            "Eg": [2, 0, -2, 0, 0, 2, 0, -2, 0, 0],
            "A1u": [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
            "A2u": [1, 1, 1, -1, -1, -1, -1, -1, 1, 1],
            "B1u": [1, -1, 1, 1, -1, -1, 1, -1, -1, 1],
            "B2u": [1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
            "Eu": [2, 0, -2, 0, 0, -2, 0, 2, 0, 0]
        }
        class_sizes = [1, 2, 1, 2, 2, 1, 2, 1, 2, 2]
        group_order = 16
        classes = ["E", "2C4", "C2", "2C2'", "2C2''", "i", "2S4", "σh", "2σv", "2σd"]
    
    elif point_group == "D5d":
        # Character table for D5d point group
        character_table = {
            "A1g": [1, 1, 1, 1, 1, 1, 1, 1],
            "A2g": [1, 1, 1, -1, 1, 1, 1, -1],
            "E1g": [2, 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), 0, 2, 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), 0],
            "E2g": [2, 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), 0, 2, 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), 0],
            "A1u": [1, 1, 1, 1, -1, -1, -1, -1],
            "A2u": [1, 1, 1, -1, -1, -1, -1, 1],
            "E1u": [2, 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), 0, -2, -2*np.cos(2*np.pi/5), -2*np.cos(4*np.pi/5), 0],
            "E2u": [2, 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), 0, -2, -2*np.cos(4*np.pi/5), -2*np.cos(2*np.pi/5), 0]
        }
        class_sizes = [1, 2, 2, 5, 1, 2, 2, 5]
        group_order = 20
        classes = ["E", "2C5", "2C5^2", "5C2", "i", "2S10", "2S10^3", "5σd"]
    
    reducible_chars = [reducible_rep[cls] for cls in classes]
    
    coefficients = {}
    for irrep, chars in character_table.items():
        coef = sum(reducible_chars[i] * chars[i] * class_sizes[i] for i in range(len(classes))) / group_order
        if np.isclose(coef, round(coef), atol=1e-10):
            coef = round(coef)
        coefficients[irrep] = coef
    
    return coefficients

def format_decomposition(coefficients):
    """Format the decomposition as a string"""
    result_str = ""
    for irrep, coef in sorted(coefficients.items()):
        if abs(coef) > 1e-10:
            if isinstance(coef, float) and np.isclose(coef, round(coef), atol=1e-10):
                coef = int(round(coef))
            if coef == 1:
                result_str += f" + {irrep}"
            elif coef > 0:
                result_str += f" + {coef}{irrep}"
            elif coef == -1:
                result_str += f" - {irrep}"
            else:
                result_str += f" - {abs(coef)}{irrep}"
    
    return result_str[3:] if result_str.startswith(" + ") else result_str or "0"

def extract_L_from_term_symbol(term_symbol):
    """Extract the L value from a term symbol (e.g., '3P' -> L=1)"""
    orbital_letter = term_symbol[1] if len(term_symbol) > 1 else term_symbol[0]
    return ORBITAL_LABELS.index(orbital_letter)

def analyze_term_symbols_in_point_group():
    """Analyze all term symbols in the user's chosen point group"""
    print("\nAnalyzing term symbols in a crystal field environment")
    print("-" * 60)
    
    # Ask user to choose point group
    print("\nChoose a point group:")
    print("1. Oh (Octahedral)")
    print("2. D4h (Tetragonal)")
    print("3. D5d (Pentagonal bipyramidal)")
    choice = input("Enter your choice (1, 2, or 3): ")
    
    if choice == "1":
        point_group = "Oh"
        calculate_reducible_rep = calculate_oh_reducible_rep
    elif choice == "2":
        point_group = "D4h"
        calculate_reducible_rep = calculate_d4h_reducible_rep
    elif choice == "3":
        point_group = "D5d"
        calculate_reducible_rep = calculate_d5d_reducible_rep
    else:
        print("Invalid choice. Defaulting to Oh point group.")
        point_group = "Oh"
        calculate_reducible_rep = calculate_oh_reducible_rep
    
    print(f"\nAnalyzing term symbols for {ORBITAL_LABELS[L]}^{r} in {point_group} point group:")
    print("-" * 60)
    
    # Extract unique L values from term symbols
    unique_L_values = set()
    for term in term_symbols_list:
        L_value = extract_L_from_term_symbol(term)
        unique_L_values.add(L_value)
    
    # Calculate and display irreducible representations for each L value
    for L_value in sorted(unique_L_values):
        orbital_label = ORBITAL_LABELS[L_value]
        print(f"\nL = {L_value} ({orbital_label} orbital):")
        
        # Calculate reducible representation
        reducible_rep = calculate_reducible_rep(L_value)
        
        # Decompose into irreducible representations
        coefficients = decompose_into_irreps(reducible_rep, point_group)
        
        # Display results
        print(f"Γ(L={L_value}) = {format_decomposition(coefficients)}")
        
        # Find which term symbols have this L value
        related_terms = [term for term in term_symbols_list if extract_L_from_term_symbol(term) == L_value]
        print(f"Related term symbols: {', '.join(related_terms)}")
        print("-" * 40)


#################################################################################
#%% CREATING HTML OUTPUT FOR WEAK FIELD ANALYSIS
def create_weak_field_html(point_group, unique_L_values, term_symbols_list):
    """Create an HTML file with the weak field analysis results"""
    # Create filename based on configuration and point group
    filename = f"weak_field_{ORBITAL_LABELS[L]}^{r}_{point_group}.html"
    
    with open(filename, 'w') as file:
        # Start HTML document
        file.write(f'''<!DOCTYPE html>
<html>
<head>
    <title>Weak Field Analysis for {ORBITAL_LABELS[L]}^{r} in {point_group}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            line-height: 1.6;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        .container {{
            max-width: 1000px;
            margin: 0 auto;
        }}
        .section {{
            margin-bottom: 30px;
            padding: 20px;
            background-color: #f9f9f9;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .term-symbol {{
            font-weight: bold;
            color: #3498db;
        }}
        .irrep {{
            font-weight: bold;
            color: #e74c3c;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .summary {{
            margin-top: 20px;
            padding: 15px;
            background-color: #e8f4f8;
            border-left: 4px solid #3498db;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Weak Field Analysis for {ORBITAL_LABELS[L]}^{r} in {point_group} Point Group</h1>
        
        <div class="section">
            <h2>Configuration Summary</h2>
            <p>Analyzing the {ORBITAL_LABELS[L]}^{r} electron configuration (L={L}, {r} electrons) in a {point_group} crystal field environment.</p>
        </div>
        
        <div class="section">
            <h2>Irreducible Representation Analysis</h2>''')
        
        # Add each L value analysis
        for L_value in sorted(unique_L_values):
            orbital_label = ORBITAL_LABELS[L_value]
            
            # Calculate reducible representation
            if point_group == "Oh":
                reducible_rep = calculate_oh_reducible_rep(L_value)
            elif point_group == "D4h":
                reducible_rep = calculate_d4h_reducible_rep(L_value)
            elif point_group == "D5d":
                reducible_rep = calculate_d5d_reducible_rep(L_value)
            
            # Decompose into irreducible representations
            coefficients = decompose_into_irreps(reducible_rep, point_group)
            
            # Find which term symbols have this L value
            related_terms = [term for term in term_symbols_list if extract_L_from_term_symbol(term) == L_value]
            
            file.write(f'''
            <div class="section">
                <h3>L = {L_value} ({orbital_label} terms)</h3>
                <p>Irreducible representation: <span class="irrep">Γ(L={L_value}) = {format_decomposition(coefficients)}</span></p>
                <p>Related term symbols: <span class="term-symbol">{', '.join(related_terms)}</span></p>
                
                <h4>Characters of the Reducible Representation</h4>
                <table>
                    <tr>
                        <th>Class</th>
                        <th>Character</th>
                    </tr>''')
            
            # Add characters table
            for cls, char in reducible_rep.items():
                if isinstance(char, float) and np.isclose(char, round(char), atol=1e-10):
                    char = int(round(char))
                file.write(f'''
                    <tr>
                        <td>{cls}</td>
                        <td>{char}</td>
                    </tr>''')
            
            file.write('''
                </table>
            </div>''')
        
        # End HTML document
        file.write('''
        </div>
        
        <div class="section">
            <h2>Interpretation</h2>
            <p>The irreducible representations above show how each atomic term symbol splits in the crystal field environment.</p>
            <p>Each term with angular momentum L splits according to its corresponding irreducible representation.</p>
            <div class="summary">
                <p>This analysis follows the weak field approximation, where the crystal field is treated as a perturbation to the free-ion states.</p>
            </div>
        </div>
    </div>
</body>
</html>''')
    
    print(f"\nWeak field analysis HTML file created: {os.path.abspath(filename)}")
    return filename

# Modify the analyze_term_symbols_in_point_group function to create the HTML output
def analyze_term_symbols_in_point_group():
    """Analyze all term symbols in the user's chosen point group and create HTML output"""
    print("\nAnalyzing term symbols in a crystal field environment")
    print("-" * 60)
    
    # Ask user to choose point group
    print("\nChoose a point group:")
    print("1. Oh (Octahedral)")
    print("2. D4h (Tetragonal)")
    print("3. D5d (Pentagonal bipyramidal)")
    choice = input("Enter your choice (1, 2, or 3): ")
    
    if choice == "1":
        point_group = "Oh"
        calculate_reducible_rep = calculate_oh_reducible_rep
    elif choice == "2":
        point_group = "D4h"
        calculate_reducible_rep = calculate_d4h_reducible_rep
    elif choice == "3":
        point_group = "D5d"
        calculate_reducible_rep = calculate_d5d_reducible_rep
    else:
        print("Invalid choice. Defaulting to Oh point group.")
        point_group = "Oh"
        calculate_reducible_rep = calculate_oh_reducible_rep
    
    print(f"\nAnalyzing term symbols for {ORBITAL_LABELS[L]}^{r} in {point_group} point group:")
    print("-" * 60)
    
    # Extract unique L values from term symbols
    unique_L_values = set()
    for term in term_symbols_list:
        L_value = extract_L_from_term_symbol(term)
        unique_L_values.add(L_value)
    
    # Calculate and display irreducible representations for each L value
    for L_value in sorted(unique_L_values):
        orbital_label = ORBITAL_LABELS[L_value]
        print(f"\nL = {L_value} ({orbital_label} orbital):")
        
        # Calculate reducible representation
        reducible_rep = calculate_reducible_rep(L_value)
        
        # Decompose into irreducible representations
        coefficients = decompose_into_irreps(reducible_rep, point_group)
        
        # Display results
        print(f"Γ(L={L_value}) = {format_decomposition(coefficients)}")
        
        # Find which term symbols have this L value
        related_terms = [term for term in term_symbols_list if extract_L_from_term_symbol(term) == L_value]
        print(f"Related term symbols: {', '.join(related_terms)}")
        print("-" * 40)
    
    # Create HTML output
    html_file = create_weak_field_html(point_group, unique_L_values, term_symbols_list)
    print(f"HTML output saved to: {html_file}")

# Run the analysis
analyze_term_symbols_in_point_group()
