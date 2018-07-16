#' ERCC ExFold RNA Spike-in Sequences
#'
#' This data contains the sequence for each of the 92 synthetic RNAs included
#' in the ERCC ExFold spike-in mix (the same RNAs are used for Mix 1 and Mix 2).
#' These spike-in mixes are the most commonly used way to included spike-in
#' control sequences into an RNA-Seq experiment. The catalog number at
#' ThermoFisher is 4456739.
#'
#' @format A \code{DNAStringSet} object (from Biostrings package) with 92 entries
#' \describe{
#'   \item{width}{the length of the sequence; accessed by
#'     Biostrings::width(ERCC92_seqs)}
#'   \item{seq}{the sequence; accessed directly by \code{as.character(ERCC92_seqs)}}
#'   \item{names}{the name of the sequence. The original names used hyphens,
#'     which are replaced with underscores ('_') here. These are the same as
#'     the rownames for \link{\code{ERCC92_data}}. Accessed by \code{names(ERCC92_seqs)}}
#'   \item{elementMetadata}{the metadata for the RNAs; accessed by
#'     \code{S4Vectors::elementMetadata(ERCC92_seq)}. This contains a data frame
#'     with 92 rows and 4 columns, with the following column names:
#'     \itemize{
#'       \item{ERCC_ID} The ID. This is the same as the names above.
#'       \item{Genbank} The Genbank entry accession number. Go to
#'         \url{https://www.ncbi.nlm.nih.gov/nuccore} to search Genbank entries.
#'       \item{X5prime_assay} The accession number for a predesigned TaqMan(R) Assay,
#'         with the amplicon in the 5' end of the synthetic RNA. Use
#'         \url{http://www.thermofisher.com/order/genome-database/browse/gene-expression/keyword/}
#'         to search for assays.
#'       \item{X3prime_assay} The accession number for a predesigned TaqMan(R) Assay,
#'         with the amplicon in the 3' end of the synthetic RNA. Use 
#'         \url{http://www.thermofisher.com/order/genome-database/browse/gene-expression/keyword/}
#'         to search for assays.
#'     }}
#' }
#' @source The product page for ERCC ExFold RNA Spike-in Mixes is found at
#'   \url{https://www.thermofisher.com/order/catalog/product/4456739}
#' @source The fasta file used to build the DNAStringSet object is found at
#'   \url{https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip}
#' @source The metadata was downloaded from
#'   \url{https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095047.txt}
"ERCC92_seqs"

#' ERCC ExFold RNA Spike-in Data
#'
#' This dataset contains the analysis data for each of the 92 synthetic RNAs
#' included in the ERCC ExFold spike-in mix (the same RNAs are used for Mix 1
#' and Mix 2). This spike-in mix is the most commonly used way to included
#' spike-in control sequences into an RNA-Seq experiment. The catalog number at
#' ThermoFisher is 4456739.
#' Each mix has four groups of RNAs with defined ratios between Mix 1 vs Mix 2:
#' Group A is 4:1; Group B is 1:1; Group C is 2:3; and Group D is 1:2. Each
#' RNA within each group has a defined concentration per microliter. The total
#' molar concentration for Mixes 1 and 2 are the same (approximately 100 nM),
#' though the total mass concentrations are slightly different.
#'
#' @format A data frame with 92 rows and 9 columns:
#' \describe{
#'   \item{rownames}{the name of the RNA (same as the sequence names in
#'     \link{\code{ERCC92_seqs}}}
#'   \item{re.sort_ID}{this is the rank order of the RNAs, sorted by group and
#'     then by most abundant}
#'   \item{subgroup}{which group is the RNA in: A, B, C, or D?}
#'   \item{length}{the length of the sequence}
#'   \item{Mix1_molar_conc}{the molar concentration of the RNA in attomoles per microliter
#'     (equivalent to picomolar molar concentration) in Mix 1. Estimated by the company.}
#'   \item{Mix1_mass_conc}{the concentration of the RNA in picograms per microliter
#'     (equivalent to micrograms per liter) in Mix 1. Calculated by the equation:
#'     mass concentration = molar concentration * 10^-6 * (transcript length *
#'     321.47 g / mol + 18.02 g / mol). The factor of 10^-6 converts the unit
#'     grams / mole to picograms / attomole to get the final unit of picograms
#'     per microliter.}
#'   \item{Mix1_TPM}{the transcripts per million value for Mix 1, calculated
#'     by the equations: copy number = molar concentration * 1 microliter *
#'     avogadro's number & TPM_i = copy number for i-th RNA / sum(copy number for all RNAs)}
#'   \item{Mix2_molar_conc}{the molar concentration of the RNA in attomoles per microliter
#'     (equivalent to picomolar molar concentration) in Mix 2. Estimated by the company.}
#'   \item{Mix2_mass_conc}{the concentration of the RNA in picograms per microliter
#'     (equivalent to micrograms per liter) in Mix 2. Calculated the same way as for Mix 1 above.}
#'   \item{Mix2_TPM}{the transcripts per million value for Mix 2, calculated the
#'     same way as for Mix 1 above.}
#'   \item{expected_foldchange_ratio}{the expected ratio between Mix 1 and Mix 2 given
#'     the subgroup of the RNA: A is 4:1, B is 1:1, C is 2:3, and D is 1:2.}
#'   \item{log2_foldchange_ratio}{the log2 transformation of the expected ratio.}
#' }
#' @source The product page for ERCC ExFold RNA Spike-in Mixes is found at
#'   \url{https://www.thermofisher.com/order/catalog/product/4456739}
#' @source The original analysis data was downloaded from
#'   \url{https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt}
#' @source The equation for calculating the molar mass of single-stranded RNA
#'   was taken from \url{https://nebiocalculator.neb.com/#!/ssrnaamt}
"ERCC92_data"
