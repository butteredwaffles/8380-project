from experiment import main
import time
from tabulate import tabulate
import random

N = 10
TRIALS = 10
ALPHABET = "abcdef"

def get_random_string():
    return ''.join(random.choices(ALPHABET, k=N)) + "$"


print(f"Performing {TRIALS} trials with {N}-length randomly generated strings with an alphabet of {ALPHABET}")

second_totals = []
for i in range(TRIALS):
    start_time = time.perf_counter()
    main(get_random_string(), print_output=False)
    finish = time.perf_counter() - start_time
    second_totals.append(finish)

print(tabulate([second_totals], showindex=True, tablefmt='double_outline'))
print("Average Seconds Taken: ", sum(second_totals) / TRIALS)



