import numpy as np

def scale_value(x, original_min, original_max, new_min, new_max):
    return new_min + ((x - original_min) * (new_max - new_min) / (original_max - original_min))


def scale_bead_chain(og_middle_points, nr_connections, new_min=1, new_max=200):
    trimmed_df = og_middle_points.head(nr_connections)
    original_min = np.min(trimmed_df)
    original_max = np.max(trimmed_df)
    df_scaled = trimmed_df.map(lambda x: scale_value(x, original_min, original_max, new_min, new_max))
    df_scaled = df_scaled.map(lambda x: int(x))
    
    return df_scaled