# O-sativa-microbiome-modeling
15 April 2020


Overall simulation and analysis workflow. 

An ensemble of ‘Native Communities’ is first generated and equilibrated. Typically this consists of initializing/randomly seeding 1000 5-member generalized Lotka Volterra systems, and finding the unique, biologically relevant fixed points. In the next step, a single ‘target species’ is introduced into the equilibrated native communities under two scenarios; inoculated alone (STSI) or  inoculated in the presence of 5 additional community members (PSC). These systems are initialized into generalized Lotka Volterra systems, and the biologically relevant fixed points are found. PSC’s which promote high abundance of the target species are further interrogated using a parameter sensitivity algorithm. 

Information about Scripts and Output Data
       
    1. Native Community Generation Script
        1. code filename: “[Date]_stode_nativeComms.R”
        2. Purpose: Generates ensemble of native communities 
        3. Returns/generates 3 .csv files:
            1. Native Community Final States – corresponding unique equilibrium states at t → inf, fixed points
                1. example filename: “NativeCommunities_stode_final_states.csv”
                2. directory that file is written to:~/Documents/SynCom_Modeling/
                3. 1000 rows, 5 columns, no headers or row labels
                4. each row - the final state (abundances) for a different community of the ensemble of native communities
                5. each of 5 columns representing one of the five total species
            2. NativeCommunity Parameters – the corresponding parameter sets which gave rise to each state, also the parameters that are carried through to subsequent simulations
                1. example filename: “NativeCommunities_stode_final_parameters.csv”
                2. directory that file is written to: ~/Documents/SynCom_Modeling/
                3. no column headers no rownames
                4. 1000 rows, 30 columns
                    1. corresponding column headers (parameter names) for this file are located in the file: …./Model_Code_Data_forAaron/ReadMes_Info/13apr20_NC_parameter_file_headers.csv  
                        1. Mu’s are growth rates (e.g. ‘muX1’ is the growth rate if species 1)
                        2. Alphas are interaction parameters (e.g. ‘a13’ is the effect species 3 has on species 1. This parameter appears as a coefficient in the equation describing the first derivative in time of the abundance of species 1 
                5. each row is a set of parameters for a different community
                6. each column a different gLV model parameter
            3. Native Community Eigenvalues - estimated at the final states for each community 
                1. example filename: “NativeCommunities_stode_final_eigs.csv”
                2. directory written to: ~/Documents/SynCom_Modeling/
                3. no column headers no rownames
                4. 1000 rows, 5 columns  
                5.  each row the set of eigenvalues for a different community
       
    2. Single Species and Probiotic Support Community Inoculation and Analysis Script
        1. code filename: “[Date]_stode_single_killerComms_analysis.R”
        2. Purpose: Generates and analyzes:
            1. Communities Resulting From Single Target Species Inoculation (STSI)
                1. Resulting in 6-species communities	
            2. Communities Resulting From Probiotic Support Communities (PSC)
                1. Resulting in 11-species communities 
        3. Each instance/run of this code generates a new directory where generated files are written to. These directories are typically numbered, and are written to the following directory: ~/Documents/SynCom_Modeling/[Directory#]/
        4. Returns/generates: 3 - directories, 8 - .csv files
            1. Directories (3): 
                1. ~/Documents/SynCom_Modeling/KillerComm_stode_community_eigs
                    1. each file is the estimated final state eigenvalues of communities resulting from an individual PSC inoculation. PSC which was inoculated is denoted in filename.
                    2. example filenames within directory: “KillerComm_[PSC_Number]_stode_final_eigs.csv”
                    3. no column headers no rownames
                    4. 1000 rows, 11 columns  
                    5. each row is the set of eigenvalues for a different resulting community
                2. ~/Documents/SynCom_Modeling/KillerComm_stode_community_parameters
                    1. each file is the set of parameters of communities (gLV systems) from an individual PSC inoculation. PSC which was inoculated is denoted in filename.
                    2. example filenames within directory: “KillerComm_[PSC_Number]_stode_final_parameters.csv”
                    3. no column headers no rownames
                    4. 1000 rows, 132 columns
                    5. each row is the set of parameters for a different resulting community
                    6. each column is a different parameter
                    7. corresponding column headers (parameter names) for this file are located in the google drive file: …./Model_Code_Data_forAaron/ReadMes_Info/13apr20_PSC_parameter_file_headers.csv  
                        1. Mu’s are growth rates (e.g. ‘muX1’ is the growth rate if species 1)
                        2. Alphas are interaction parameters (e.g. ‘a13’ is the effect species 3 has on species 1. This parameter appears as a coefficient in the equation describing the first derivative in time of the abundance of species 1 
                3. ~/Documents/SynCom_Modeling/KillerComm_stode_community_states
                    1. each file is the set of states of communities resulting from an individual PSC inoculation. PSC which was inoculated is denoted in filename
                    2. example filenames within directory: “KillerComm_[PSC_Number]_stode_final_states.csv”
                    3. each file has no column headers no rownames
                    4. each file has 1000 rows, 11 columns  
                    5. each row is the set of final states for a different resulting community
                    6. each column represents a different species/states
            2. .csv files (8):
                1. Native Community Final States
                    1. Copied from files generated in step 1 (Native Community Generation), denoting which native communities were used in STSI and PSC inoculation
                2. Native Community Parameters 
                    1. Copied from files generated in step 1 (Native Community Generation), denoting which native communities were used in STSI and PSC inoculation
                3. Native Community Eigenvalueshe above 3 files are the same files 
                    1. Copied from files generated in step 1 (Native Community Generation), denoting which native communities were used in STSI and PSC inoculation
                4. Single Target Species Inoculation Final States
                    1. estimated final state eigenvalues of communities resulting from Single Target Species Inoculation (STSI).
                    2. example filename: “SingleKiller_stode_final_states.csv”
                    3. each file has no column headers no rownames
                    4. each file has 1000 rows, 6 columns  
                    5. each row is the set of final states for a different resulting community
                    6. each column represents a different individual species
                5. Single Target Species Inoculation Parameters
                    1. each file is the set of parameters of communities resulting from an individual PSC inoculation. PSC which was inoculated is denoted in filename
                    2. example filenames within directory: “SingleKiller_stode_final_parameters.csv”
                    3. no column headers no rownames
                    4. 1000 rows, 42 columns  
                    5. each row is the set of parameters for a different resulting community
                    6. each column is a different parameter
                    7. corresponding column headers (parameter names) for this file are located in the file: …./Model_Code_Data_forAaron/ReadMes_Info/13apr20_STSI_parameter_file_headers.csv  
                        1. Mu’s are growth rates (e.g. ‘muX1’ is the growth rate if species 1)
                        2. Alphas are interaction parameters (e.g. ‘a13’ is the effect species 3 has on species 1. This parameter appears as a coefficient in the equation describing the first derivative in time of the abundance of species 1 
                6. Single Target Species Inoculation Eigenvalues
                    1. each file is the estimated final state eigenvalues of communities resulting from an individual STSI inoculation. 
                    2. example filenames within directory: “SingleKiller_stode_final_eigs.csv”
                    3. no column headers no rownames
                    4. 1000 rows, 6 columns  
                    5. each row is the set of eigenvalues for a different resulting community
                7. Combined PSC Metrics
                    1. Example Filename: CombinedMetrics_file_number_[file#].csv
                    2. Aggregate metrics for whole set of PSC’s generated in particular set
                    3. Single row of data, with column header labels
                8. Individual PSC Metrics
                    1. Example Filename: KillerCommMetrics_file_number_[file#].csv
                    2. Aggregate metrics for each PSC, across all resulting communities generated for each PSC inoculation
                    3. Number of rows is number of PSC’s generated and surveyed (usually between 200 and 2000)
                    4. Columns are metrics/data which are computed for each PSC, columns have header labels
       
    3. Parameter Sensitivity Analysis – Functional Motif Identification Script
        1. code filename: “[Date]_parameterSensitivity_motifFinding.R”
        2. Purpose: Analyzes individual parameter sensitivity for a chosen/given PSC. Does so by, individually and separately, setting each parameter to zero (removes graph weight/edge) and assessing modified community metrics with that parameter omitted  
        3. User provides PSC community number (and corresponding set number) to analyze in code file header
        4. Returns/generates - 1 x .csv files
            1. Parameter Sensitivity Metrics
                1. Example filename: ParameterSensitivity_PSC_[PSC#].csv
                2. 80 rows – each row is a ‘common’ parameter being inividually omitted
                3. Columns are metrics being computed for each parameter omission. Columns have header labels.
