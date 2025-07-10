
def format_value_ci(value, ci_lower, ci_upper, sig_digits=2):
    """
    Format value and 95% CI bounds per Hughes and Hase significant digit rules.
    
    Parameters:
    - value: Measured value (e.g., proportion 0.783648374)
    - ci_lower, ci_upper: 95% CI bounds
    - sig_digits: Significant digits for uncertainty (1 or 2)
    
    Returns:
    - Formatted string with value, uncertainty, percentage uncertainty, and CI
    """
    # Compute uncertainty (half CI width)
    uncertainty = (ci_upper - ci_lower) / 2
    
    # Round uncertainty to sig_digits
    if sig_digits == 1:
        uncertainty = round(uncertainty, -int(np.floor(np.log10(uncertainty))) + 1)
    else:
        uncertainty = round(uncertainty, -int(np.floor(np.log10(uncertainty))) + 2)
    
    # Determine decimal places from uncertainty
    decimal_places = -int(np.floor(np.log10(uncertainty)))
    
    # Round value and CI bounds
    value = round(value, decimal_places)
    ci_lower = round(ci_lower, decimal_places)
    ci_upper = round(ci_upper, decimal_places)
    
    # Percentage uncertainty
    rel_uncertainty = uncertainty / value
    pct_uncertainty = round(rel_uncertainty * 100, 1)  # Two sig digits
    
    # Format output
    fmt = f".{decimal_places}f"
    return (f"Value: {value:{fmt}} ± {uncertainty:{fmt}} "
            f"(95% CI: [{ci_lower:{fmt}}, {ci_upper:{fmt}}])\n"
            f"Percentage: {value*100:{fmt.replace('f', '0f')}}% ± {pct_uncertainty:.1f}% "
            f"(95% CI: [{ci_lower*100:{fmt.replace('f', '0f')}}%, {ci_upper*100:{fmt.replace('f', '0f')}}%])")


if __name__ == "__main__":
    # Example usage
    import numpy as np

    # Example values from Hughes and Hase (2019)
    # These values are illustrative; replace with actual data as needed.
    # Value, CI lower, CI upper

    # Your data
    value = 0.02644845324359870048
    ci_lower, ci_upper = 0.00000000015707745481, 0.05324877359942899874
    
    # Format and print the value with CI
    value = 0.73437110994498255856
    ci_lower = 0.02246928130054030165
    ci_upper = 0.91242909875254762930

    print(format_value_ci(value, ci_lower, ci_upper, sig_digits=2))