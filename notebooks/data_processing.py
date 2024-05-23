import numpy as np

def scale_bead_chain(og_middle_points, nr_connections, new_min=1, new_max=200):
    trimmed_df = og_middle_points.head(nr_connections)
    original_min = np.min(trimmed_df)
    original_max = np.max(trimmed_df)
    df_scaled = new_min + ((trimmed_df - original_min) * (new_max - new_min) / (original_max - original_min))
    
    return df_scaled.astype('int64')