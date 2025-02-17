import csv
from collections import defaultdict
from itertools import combinations

# Initialize data structures
category_counts = defaultdict(set)
rows_by_category = defaultdict(list)

# Read the input file
with open('output_final_NEW.txt', 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        id_ = row[0]
        category = row[1]
        category_counts[id_].add(category)
        rows_by_category[id_].append(row)

# Define all possible categories
all_categories = {'P1', 'P3', 'Tpase', 'Mutator'}

# Initialize counters for statistics and files for each combination
combination_counts = defaultdict(int)
files = {}

# Generate file handles for all combinations
for i in range(1, len(all_categories) + 1):
    for combo in combinations(all_categories, i):
        combo_name = '_'.join(sorted(combo))
        files[combo_name] = open(f'{combo_name}.txt', 'w')

# Calculate statistics and write rows to corresponding files
total_ids = len(category_counts)

for id_, categories in category_counts.items():
    sorted_categories = sorted(categories)
    combo_name = '_'.join(sorted_categories)
    combination_counts[combo_name] += 1
    for row in rows_by_category[id_]:
        files[combo_name].write('\t'.join(row) + '\n')
# Close all files
for f in files.values():
    f.close()

# Print results
print(f'Total IDs: {total_ids}')
for combo, count in combination_counts.items():
    percentage = (count / total_ids) * 100
    print(f'{combo}: {count} ({percentage:.2f}%)')
