import numpy as np

def load_csv_data(in_file: str) -> None:
    data = np.loadtxt(in_file, delimiter=',', dtype=float)
    return data


if __name__ == '__main__':
    print(load_csv_data("efficiency_data/compound_data/CsI.csv"))  