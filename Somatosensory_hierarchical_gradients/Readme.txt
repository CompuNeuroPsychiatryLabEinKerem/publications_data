These codes are part of our publication, "Hierarchical cortical gradients in somatosensory processing". 
https://doi.org/10.1016/j.neuroimage.2020.117257
For questions please contact Noam Saadon-Grosman, noam.saadongros@mail.huji.ac.il.
We used BrainVoyager, Matlab  and Neuroelf for analysis. Other formats can be used as Matlab matrices, in that case- the xff command should be excluded. 

Hierarchical_gradients_1_Somatosensory_response_tMap- takes mesh (surface) time courses of all participants and the two stimulus directions (star_lip and start_toe) and computes a random effect group t-values map corrected for multiple comparisons and masked for significance.

Hierarchical_gradients_2_Selectivity_computation- computes selectivity in each vertex of each subject in each stimulus direction that passed the threshold for significant somatosensory response. Selectivity is defined as 1/standard deviation of a Gaussian fitted to the response Event Related Averaging.

Hierarchical_gradients_3_Laterality_computation- computes laterality in each vertex of each subject in each stimulus direction that passed the threshold for significant somatosensory response. Laterality is defined as the normalized difference between contralateral and ipsilateral response.

Hierarchical_gradients_4_Selectivity_map- creates selectivity surface map from the output of "Selectivity_computation"

Hierarchical_gradients_5_Laterality_map- creates laterality surface map from the output of "Laterality_computation"

Hierarchical_gradients_6_geodesic_distance_to_central_sulcus- This code computes the geodesic distance between each vertex to the closest vertex in area 3a. We use an online code (https://code.google.com/archive/p/geodesic/) to compute the distance, separately for each gross anatomical region. It also computes selectivity and laterality values within distance quantiles.
