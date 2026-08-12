# Applied in order to title-cased hallmark names: acronyms first, then words
PATHWAY_ABBREV <- c(
  "Epithelial Mesenchymal Transition" = "EMT",
  "Reactive Oxygen Species"           = "ROS",
  "Unfolded Protein Response"         = "UPR",
  "Oxidative Phosphorylation"         = "OxPhos",
  "Tnfa"                              = "TNFA",
  "Nfkb"                              = "NF-kB",
  "Il6"                               = "IL6",
  "Il2"                               = "IL2",
  "Jak"                               = "JAK",
  "Stat3"                             = "STAT3",
  "Stat5"                             = "STAT5",
  "E2f"                               = "E2F",
  "G2m"                               = "G2M",
  "Myc"                               = "MYC",
  "Kras"                              = "KRAS",
  "Mtorc1"                            = "mTORC1",
  "Pi3k"                              = "PI3K",
  "Akt"                               = "AKT",
  "Tgf"                               = "TGF",
  "Wnt"                               = "WNT",
  "Dna"                               = "DNA",
  "Uv"                                = "UV",
  "P53"                               = "p53",
  "Interferon"                        = "IFN",
  "Signaling"                         = "Sig.",
  "Response"                          = "Resp.",
  "Metabolism"                        = "Metab.",
  "Via"                               = "via"
)

# Helper to shorten hallmark pathway names for axis labels
abbrev_pathway <- function(x, max_chars = 30) {
  x <- str_to_title(str_replace_all(str_remove(x, "^HALLMARK_"), "_", " "))
  x <- str_replace_all(x, PATHWAY_ABBREV)
  str_trunc(x, max_chars)
}
