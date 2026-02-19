import numpy as np

def save_columns_to_file(
    output_data_filepath: str, 
    column_width: int, 
    precision: int, 
    tuple_of_columns: list,
    column_labels: list[str]
) -> None:
    
    data = np.column_stack(tuple_of_columns)
    
    column_labels_with_spaces = []
    for label in column_labels:
        empty_spaces = column_width - len(label)
        for _ in range(empty_spaces):
            label = " " + label
        column_labels_with_spaces.append(label)

    header = " ".join(column_labels_with_spaces)

    np.savetxt(
        output_data_filepath,
        data,
        header=header,
        comments="",
        fmt=f'%{column_width}.{precision}f',
    )
