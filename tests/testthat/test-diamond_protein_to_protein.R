context("Test: diamond_protein_to_protein()")

test_that("The diamond example works properly", {

  skip_on_cran()
  skip_on_travis()

  diamond_example <- diamond_protein_to_protein(
    query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
    subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
    sensitivity_mode = "ultra-sensitive",
    diamond_exec_path = "/opt/miniconda3/bin/",
    output_path = tempdir(),
    db_import  = FALSE,
    cores = 1
  )

  expect_equal(nrow(diamond_example), 20)
  expect_equal(ncol(diamond_example), 20)


})

test_that("All diamond sensitivity modes work properly", {
  skip_on_cran()
  skip_on_travis()

  # test mode: fast
    diamond_example_fast <- diamond_protein_to_protein(
      query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
      subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
      sensitivity_mode = "fast",
      diamond_exec_path = "/opt/miniconda3/bin/",
      output_path = tempdir(),
      db_import  = FALSE,
      cores = 1
    )
    expect_equal(nrow(diamond_example_fast), 20)
    expect_equal(ncol(diamond_example_fast), 20)

    # test mode: mid-sensitive
    diamond_example_mid_sensitive <- diamond_protein_to_protein(
      query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
      subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
      sensitivity_mode = "mid-sensitive",
      diamond_exec_path = "/opt/miniconda3/bin/",
      output_path = tempdir(),
      db_import  = FALSE,
      cores = 1
    )

    expect_equal(nrow(diamond_example_mid_sensitive), 20)
    expect_equal(ncol(diamond_example_mid_sensitive), 20)

    # test mode: sensitive
    diamond_example_sensitive <- diamond_protein_to_protein(
      query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
      subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
      sensitivity_mode = "sensitive",
      diamond_exec_path = "/opt/miniconda3/bin/",
      output_path = tempdir(),
      db_import  = FALSE,
      cores = 1
    )

    expect_equal(nrow(diamond_example_sensitive), 20)
    expect_equal(ncol(diamond_example_sensitive), 20)

  # test mode: more-sensitive
    diamond_example_more_sensitive <- diamond_protein_to_protein(
      query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
      subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
      sensitivity_mode = "more-sensitive",
      diamond_exec_path = "/opt/miniconda3/bin/",
      output_path = tempdir(),
      db_import  = FALSE,
      cores = 1
    )

    expect_equal(nrow(diamond_example_more_sensitive), 20)
    expect_equal(ncol(diamond_example_more_sensitive), 20)


  # test mode: very-sensitive
    diamond_example_very_sensitive <- diamond_protein_to_protein(
      query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
      subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
      sensitivity_mode = "very-sensitive",
      diamond_exec_path = "/opt/miniconda3/bin/",
      output_path = tempdir(),
      db_import  = FALSE,
      cores = 1
    )

    expect_equal(nrow(diamond_example_very_sensitive), 20)
    expect_equal(ncol(diamond_example_very_sensitive), 20)


  # test mode: ultra-sensitive
    diamond_example_ultra_sensitive <- diamond_protein_to_protein(
      query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
      subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
      sensitivity_mode = "ultra-sensitive",
      diamond_exec_path = "/opt/miniconda3/bin/",
      output_path = tempdir(),
      db_import  = FALSE,
      cores = 1
    )

    expect_equal(nrow(diamond_example_ultra_sensitive), 20)
    expect_equal(ncol(diamond_example_ultra_sensitive), 20)
})

test_that("diamond runs properly with various max-target-seqs options", {
  skip_on_cran()
  skip_on_travis()

  diamond_example_ultra_sensitive <- diamond_protein_to_protein(
    query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
    subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
    sensitivity_mode = "ultra-sensitive",
    diamond_exec_path = "/opt/miniconda3/bin/",
    max_target_seqs = 1,
    output_path = tempdir(),
    db_import  = FALSE,
    cores = 1
  )

  expect_equal(nrow(diamond_example_ultra_sensitive), 20)
  expect_equal(ncol(diamond_example_ultra_sensitive), 20)
})



test_that("diamond runs properly when additional parameter options are passed along", {
  skip_on_cran()
  skip_on_travis()

diamond_example_ultra_sensitive_add_diamond_options <- diamond_protein_to_protein(
  query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
  subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
  sensitivity_mode = "ultra-sensitive",
  diamond_exec_path = "/opt/miniconda3/bin/",
  max_target_seqs = 1,
  output_path = tempdir(),
  db_import  = FALSE,
  add_diamond_options = "--block-size 4.0 --compress 1 --no-self-hits",
  cores = 1
)

expect_equal(nrow(diamond_example_ultra_sensitive_add_diamond_options), 20)
expect_equal(ncol(diamond_example_ultra_sensitive_add_diamond_options), 20)


diamond_example_ultra_sensitive_add_makedb_options <- diamond_protein_to_protein(
  query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
  subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
  sensitivity_mode = "ultra-sensitive",
  diamond_exec_path = "/opt/miniconda3/bin/",
  max_target_seqs = 1,
  output_path = tempdir(),
  db_import  = FALSE,
  add_makedb_options = "--taxonnames",
  cores = 1
)

expect_equal(nrow(diamond_example_ultra_sensitive_add_makedb_options), 20)
expect_equal(ncol(diamond_example_ultra_sensitive_add_makedb_options), 20)


})
