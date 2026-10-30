CONFIG <- list(
  project_name = "PRJNA932038_ASC_P_human_validation",
  input_root = "D:/Users/lstsw/Documents/ASC_CD4/final_result/final_result",
  seed = 1234L,
  samples = data.frame(
    folder = c("s1_SRR23526669", "s2_SRR23526668", "s3_SRR23526667",
               "s4_SRR23526666", "s5_SRR23526665", "s6_SRR23526664"),
    run_accession = c("SRR23526669", "SRR23526668", "SRR23526667",
                      "SRR23526666", "SRR23526665", "SRR23526664"),
    patient_id = c("S53", "S53", "S76", "S76", "S81", "S81"),
    group = c("N", "T", "N", "T", "N", "T"),
    sample_id = c("S53_N", "S53_T", "S76_N", "S76_T", "S81_N", "S81_T"),
    stringsAsFactors = FALSE
  ),
  qc = list(
    min_features = 200L,
    max_percent_mt = 10,
    mad_multiplier = 3,
    article_sd_multiplier = 2
  ),
  integration = list(
    npcs = 30L,
    dims = 1:30,
    global_resolution = 0.4,
    stromal_resolution = 0.4,
    harmony_group = "sample_id"
  ),
  signature = c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM"),
  tcell_gate = list(min_total = 50L, min_per_patient = 10L, min_patients = 2L),
  pseudobulk = list(min_cells_per_sample_status = 20L, fdr = 0.05),
  local_library = "R_library"
)

