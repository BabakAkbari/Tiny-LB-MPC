import numpy as np
import pandas as pd

def format_matrix(name, matrix):
    """Format a NumPy matrix as an Eigen::Matrix declaration."""
    rows, cols = matrix.shape
    matrix_str = ",\n".join(
        ",".join(f"{val:.18e}" for val in row) for row in matrix
    )
    return f"static Eigen::Matrix<double, {rows}, {cols}> {name} = (Eigen::Matrix<double, {rows}, {cols}>() <<\n{matrix_str}\n).finished();\n"

def generate_gp_h(x_file, y_file, output_file="gp.h"):
    """Generate gp.h from X_train_sampled.csv and Y_train_sampled.csv."""
    X = pd.read_csv(x_file, header=None).values
    Y = pd.read_csv(y_file, header=None).values
    
    with open(output_file, "w") as f:
        f.write("#include <Eigen/Dense>\n\n")
        f.write(format_matrix("X_x", X))
        f.write("\n")
        f.write(format_matrix("Y_x", Y))
        f.write("\n")
        f.write(format_matrix("X_y", X))
        f.write("\n")
        f.write(format_matrix("Y_y", Y))
        f.write("\n")

if __name__ == "__main__":
    generate_gp_h("X_train_sampled.csv", "Y_train_sampled.csv")
