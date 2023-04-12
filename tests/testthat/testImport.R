test_that("import and evaluate DB", {
  db <- import_db("exampleFasta.fasta")
  expect_equal(attr(db[[1]], "name"), "tr|A0A075B6G3|A0A075B6G3_HUMAN")
  expect_equal(length(db), 75)
  expect_equal(nchar(db[[5]]), 551)
  protein <- select_prot(db, "A0A1B0GWI1")
  expect_error(select_prot(db, "AAAAAAAAAAAA"))
  expect_equal(protein[[1]], "MTSGANSSGSYLPSEIRSSKIDDNYLKELNEDLKLRKQELLEMLKPLEDKNNLLFQKLMSNLEEKQRSLQIMRQIMAGKGCEESSVMELLKEAEEMKQNLERKNKMLRKEMEMLWNKTFEAEELSDQQKAPQTKNKADLQDGKEKQQRKMEWVKYQEQNNILQNDFHGKVIELRIEALKNYQKANDLKLSLYLQQNFEPMQAFLNLPGSQGIRIV")
  protein_volcano <- select_prot_volcano(db, "A0A1B0GWI1")
  expect_equal(protein_volcano, "A0A1B0GWI1")
  protein_attributes <- protein_attributes(protein)
  expect_equal(protein_attributes, c("MTSGANSSGSYLPSEIRSSKIDDNYLKELNEDLKLRKQELLEMLKPLEDKNNLLFQKLMSNLEEKQRSLQIMRQIMAGKGCEESSVMELLKEAEEMKQNLERKNKMLRKEMEMLWNKTFEAEELSDQQKAPQTKNKADLQDGKEKQQRKMEWVKYQEQNNILQNDFHGKVIELRIEALKNYQKANDLKLSLYLQQNFEPMQAFLNLPGSQGIRIV",
                                     "tr|A0A1B0GWI1|A0A1B0GWI1_HUMAN",
                                     ">tr|A0A1B0GWI1|A0A1B0GWI1_HUMAN Putative coiled-coil domain-containing protein 196 OS=Homo sapiens OX=9606 GN=CCDC196 PE=1 SV=1",
                                     "A0A1B0GWI1"))
  AA_df <- create_AA_df(protein)
  expect_true(is.data.frame(AA_df))
  expect_equal(AA_df[[1]][45], "K")

})


test_that("import PEAKS", {
  PEAKS_peptides_separate <- read_peptide_csv_PEAKS_bysamp("PEAKS.exampleBySamp.peptides.csv")
  expect_equal(PEAKS_peptides_separate[[1]]$Peptide[4], "AALTELSLGSAYQAM(+15.99)ILGVDSK")
  expect_equal(nrow(PEAKS_peptides_separate[[1]]), 13)
  expect_equal(sum(grepl("PSM", names(PEAKS_peptides_separate[[1]]))), 1)
  expect_equal(sum(PEAKS_peptides_separate[[1]]$Decoy), 0)

  #PEAKS_peptides_comb <- read_peptide_csv_PEAKS_comb()
})
