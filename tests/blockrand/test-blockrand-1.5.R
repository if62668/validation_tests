# Tests for package blockrand version 1.5

library(blockrand)
library(randomizr)
library(testthat)


# Basic two-arm randomisation with actual block size 2
# With two treatment arms, block.sizes = 1 corresponds
# to an actual block size of 2.

set.seed(123)

x2 <- blockrand(
  n = 20,
  levels = c("A", "B"),
  block.sizes = 1
)

test_that("blockrand returns the expected output structure", {
  expect_s3_class(x2, "data.frame")
  
  expect_true(
    all(c("id", "block.id", "block.size", "treatment") %in% names(x2))
  )
  
  expect_equal(nrow(x2), 20)
  expect_equal(length(unique(x2$id)), nrow(x2))
})

test_that("actual block size 2 gives one allocation to each treatment", {
  blocks <- split(x2, x2$block.id)
  
  for (b in blocks) {
    expect_equal(
      nrow(b),
      2
    )
    
    expect_equal(
      unique(b$block.size),
      2
    )
    
    tab <- table(b$treatment)
    
    expect_equal(
      as.integer(tab[c("A", "B")]),
      c(1L, 1L)
    )
  }
})


# Basic two-arm randomisation with actual block size 4
# With two treatment arms, block.sizes = 2 corresponds
# to an actual block size of 4.

set.seed(123)

x4 <- blockrand(
  n = 20,
  levels = c("A", "B"),
  block.sizes = 2
)

test_that("actual block size 4 gives two allocations to each treatment", {
  blocks <- split(x4, x4$block.id)
  
  for (b in blocks) {
    expect_equal(
      nrow(b),
      4
    )
    
    expect_equal(
      unique(b$block.size),
      4
    )
    
    tab <- table(b$treatment)
    
    expect_equal(
      as.integer(tab[c("A", "B")]),
      c(2L, 2L)
    )
  }
})


# Varying actual block sizes of 2 and 4
# With two treatment arms, block.sizes = 1:2 corresponds
# to actual block sizes of 2 and 4.

set.seed(123)

x_var <- blockrand(
  n = 120,
  levels = c("A", "B"),
  block.sizes = 1:2
)

test_that("only actual block sizes 2 and 4 are generated", {
  expect_true(
    all(x_var$block.size %in% c(2, 4))
  )
})

test_that("both actual block sizes 2 and 4 occur with the fixed seed", {
  expect_equal(
    sort(unique(x_var$block.size)),
    c(2, 4)
  )
})

test_that("actual block size equals the number of allocations in each block", {
  blocks <- split(x_var, x_var$block.id)
  
  for (b in blocks) {
    expect_length(
      unique(b$block.size),
      1
    )
    
    expect_equal(
      nrow(b),
      unique(b$block.size)
    )
  }
})

test_that("varying block sizes 2 and 4 maintain 1:1 allocation within each block", {
  blocks <- split(x_var, x_var$block.id)
  
  for (b in blocks) {
    tab <- table(b$treatment)
    
    expect_equal(
      as.integer(tab[c("A", "B")]),
      rep(nrow(b) / 2, 2)
    )
  }
})


# Additional systematic check of block.sizes = 1, 2 and 3
# For two treatment arms, these correspond to actual
# block sizes 2, 4 and 6.

test_that("block.sizes 1, 2 and 3 produce correctly balanced blocks", {
  for (bs in 1:3) {
    set.seed(123)
    
    x <- blockrand(
      n = 60,
      levels = c("A", "B"),
      block.sizes = bs
    )
    
    blocks <- split(x, x$block.id)
    expected_block_size <- 2 * bs
    
    for (b in blocks) {
      expect_equal(
        nrow(b),
        expected_block_size
      )
      
      expect_equal(
        unique(b$block.size),
        expected_block_size
      )
      
      tab <- table(b$treatment)
      
      expect_equal(
        as.integer(tab[c("A", "B")]),
        c(bs, bs)
      )
    }
  }
})


# Stratified block randomisation
#
# Four strata are randomised independently.
# Within each stratum, two treatment arms and varying
# actual block sizes of 2 and 4 are used.
#
# blockrand completes the final block. Therefore, requesting
# n = 120 does not necessarily result in exactly 120 positions
# in the generated randomisation list.

set.seed(123)

total.n <- 120
groups <- c("A", "B")

stratumA <- blockrand(
  n = total.n,
  num.levels = length(groups),
  levels = groups,
  block.sizes = 1:2,
  stratum = "S1",
  id.prefix = "A",
  block.prefix = "A"
)

stratumB <- blockrand(
  n = total.n,
  num.levels = length(groups),
  levels = groups,
  block.sizes = 1:2,
  stratum = "S2",
  id.prefix = "B",
  block.prefix = "B"
)

stratumC <- blockrand(
  n = total.n,
  num.levels = length(groups),
  levels = groups,
  block.sizes = 1:2,
  stratum = "S3",
  id.prefix = "C",
  block.prefix = "C"
)

stratumD <- blockrand(
  n = total.n,
  num.levels = length(groups),
  levels = groups,
  block.sizes = 1:2,
  stratum = "S4",
  id.prefix = "D",
  block.prefix = "D"
)

rand <- rbind(
  stratumA,
  stratumB,
  stratumC,
  stratumD
)


# Check that all four strata are generated.

test_that("all four strata are generated", {
  expect_equal(
    sort(unique(as.character(rand$stratum))),
    c("S1", "S2", "S3", "S4")
  )
})


# blockrand completes the final block, so each generated
# randomisation list must contain at least 120 positions.

test_that("generated randomisation lists contain at least 120 positions per stratum", {
  n_per_stratum <- table(rand$stratum)
  
  expect_true(
    all(n_per_stratum >= total.n)
  )
})


# Check actual block sizes across all four strata.

test_that("only actual block sizes 2 and 4 occur across strata", {
  expect_true(
    all(rand$block.size %in% c(2, 4))
  )
})

test_that("both actual block sizes 2 and 4 occur across strata with the fixed seed", {
  expect_equal(
    sort(unique(rand$block.size)),
    c(2, 4)
  )
})


# Check internal consistency of every block in every stratum.

test_that("each block is internally consistent within every stratum", {
  blocksA <- split(stratumA, stratumA$block.id)
  blocksB <- split(stratumB, stratumB$block.id)
  blocksC <- split(stratumC, stratumC$block.id)
  blocksD <- split(stratumD, stratumD$block.id)
  
  for (b in c(blocksA, blocksB, blocksC, blocksD)) {
    expect_length(
      unique(b$block.size),
      1
    )
    
    expect_equal(
      nrow(b),
      unique(b$block.size)
    )
    
    expect_true(
      unique(b$block.size) %in% c(2, 4)
    )
  }
})


# Check 1:1 treatment allocation within every complete block.
#
# Actual block size 2: one allocation to A and one to B.
# Actual block size 4: two allocations to A and two to B.

test_that("treatment allocation is 1:1 within every complete block and stratum", {
  blocksA <- split(stratumA, stratumA$block.id)
  blocksB <- split(stratumB, stratumB$block.id)
  blocksC <- split(stratumC, stratumC$block.id)
  blocksD <- split(stratumD, stratumD$block.id)
  
  for (b in c(blocksA, blocksB, blocksC, blocksD)) {
    tab <- table(b$treatment)
    
    expect_equal(
      as.integer(tab[groups]),
      rep(nrow(b) / 2, 2)
    )
  }
})


# Allocation at the requested sample size
#
# blockrand completes the final block and may therefore generate
# more positions than requested. Only the first 120 positions
# within each stratum are considered here.

usedA <- stratumA[seq_len(total.n), ]
usedB <- stratumB[seq_len(total.n), ]
usedC <- stratumC[seq_len(total.n), ]
usedD <- stratumD[seq_len(total.n), ]

used <- rbind(
  usedA,
  usedB,
  usedC,
  usedD
)

test_that("exactly 120 allocations per stratum are used", {
  expect_true(
    all(table(used$stratum) == total.n)
  )
})


# Complete blocks are exactly balanced.
# If the first 120 positions cut through the final block,
# a temporary treatment imbalance can occur.
# With a largest actual block size of 4, the treatment-count
# difference is checked to be no greater than 2.

test_that("treatment imbalance at n = 120 is limited within each stratum", {
  tab <- table(
    used$stratum,
    used$treatment
  )
  
  imbalance <- abs(
    tab[, groups[1]] - tab[, groups[2]]
  )
  
  expect_true(
    all(imbalance <= 2)
  )
})


# Reproducibility
#
# Identical seeds and identical arguments should produce
# identical randomisation lists.


# Reproducibility for actual block size 2.

test_that("randomisation with block size 2 is reproducible with the same seed", {
  set.seed(123)
  
  x1 <- blockrand(
    n = 60,
    levels = c("A", "B"),
    block.sizes = 1
  )
  
  set.seed(123)
  
  x2 <- blockrand(
    n = 60,
    levels = c("A", "B"),
    block.sizes = 1
  )
  
  expect_identical(x1, x2)
})


# Reproducibility for actual block size 4.

test_that("randomisation with block size 4 is reproducible with the same seed", {
  set.seed(123)
  
  x1 <- blockrand(
    n = 60,
    levels = c("A", "B"),
    block.sizes = 2
  )
  
  set.seed(123)
  
  x2 <- blockrand(
    n = 60,
    levels = c("A", "B"),
    block.sizes = 2
  )
  
  expect_identical(x1, x2)
})


# Reproducibility for varying actual block sizes 2 and 4.

test_that("randomisation with varying block sizes 2 and 4 is reproducible with the same seed", {
  set.seed(123)
  
  x1 <- blockrand(
    n = 120,
    levels = c("A", "B"),
    block.sizes = 1:2
  )
  
  set.seed(123)
  
  x2 <- blockrand(
    n = 120,
    levels = c("A", "B"),
    block.sizes = 1:2
  )
  
  expect_identical(x1, x2)
})


# Independent comparison with randomizr
#
# The same design requirements are specified independently
# for blockrand and randomizr.
#
# The exact random treatment sequence is not compared because
# the two packages use independent randomisation procedures.
# Instead, the treatment counts within each complete block
# are compared.


# Independent comparison for actual block size 2.

test_that("block size 2 agrees with independent implementation", {
  set.seed(123)
  
  br <- blockrand(
    n = 20,
    levels = c("A", "B"),
    block.sizes = 1
  )
  
  rz_blocks <- rep(seq_len(10), each = 2)
  
  set.seed(456)
  
  rz <- block_ra(
    blocks = rz_blocks,
    block_m_each = matrix(
      rep(c(1, 1), 10),
      ncol = 2,
      byrow = TRUE
    ),
    conditions = c("A", "B")
  )
  
  br_counts <- table(
    br$block.id,
    br$treatment
  )
  
  rz_counts <- table(
    rz_blocks,
    rz
  )
  
  expect_equal(
    unname(br_counts),
    unname(rz_counts)
  )
})


# Independent comparison for actual block size 4.

test_that("block size 4 agrees with independent implementation", {
  set.seed(123)
  
  br <- blockrand(
    n = 20,
    levels = c("A", "B"),
    block.sizes = 2
  )
  
  rz_blocks <- rep(seq_len(5), each = 4)
  
  set.seed(456)
  
  rz <- block_ra(
    blocks = rz_blocks,
    block_m_each = matrix(
      rep(c(2, 2), 5),
      ncol = 2,
      byrow = TRUE
    ),
    conditions = c("A", "B")
  )
  
  br_counts <- table(
    br$block.id,
    br$treatment
  )
  
  rz_counts <- table(
    rz_blocks,
    rz
  )
  
  expect_equal(
    unname(br_counts),
    unname(rz_counts)
  )
})