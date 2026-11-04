# ============================================================================ #
# test-00-iwmde-row-state-groups.R
# ============================================================================ #

context("IWMDE active-branch grouping cache")
skip_on_cran()


# A minimal stand-in for the row states the grouping reads: only the row index
# and the field .iwmde_state_active_key() resolves are needed here.
.group_state <- function(row_index, key) {

  return(list(row_index = row_index, active_key = key))
}


.group_context <- function() {

  return(list(active_key_stub = TRUE))
}


test_that("the cached grouping is rejected when a row index changes", {

  local_mocked_bindings(
    .iwmde_state_active_key = function(context, state) state[["active_key"]],
    .package = "RoBMA"
  )

  states <- list(
    .group_state(3L, "a"),
    .group_state(7L, "b"),
    .group_state(11L, "a")
  )
  context <- .group_context()

  named <- function(x) {
    attributes(x) <- list(names = names(x))
    x
  }

  groups <- .iwmde_row_state_groups(context, states)
  expect_identical(named(groups), list(a = c(1L, 3L), b = 2L))
  expect_identical(attr(groups, "n_states", exact = TRUE), 3L)
  expect_identical(attr(groups, "row_edges", exact = TRUE), c(3L, 11L))

  # Serving the cache when the states are unchanged.
  attr(states, "iwmde_active_groups") <- groups
  expect_identical(.iwmde_row_state_groups(context, states), groups)

  # A stale grouping whose length still matches is rejected at either edge.
  edited <- states
  edited[[1L]] <- .group_state(4L, "c")
  expect_identical(named(.iwmde_row_state_groups(context, edited)),
                   list(c = 1L, b = 2L, a = 3L))

  edited <- states
  edited[[3L]] <- .group_state(12L, "a")
  refreshed <- .iwmde_row_state_groups(context, edited)
  expect_identical(attr(refreshed, "row_edges", exact = TRUE), c(3L, 12L))

  # An interior edit that keeps both edges still yields the cached grouping;
  # the guard covers in-place length-preserving edits at the edges, which is
  # what the subsetting call sites can produce.
  edited <- states
  edited[[2L]] <- .group_state(7L, "b")
  expect_identical(.iwmde_row_state_groups(context, edited), groups)

  # Subsetting drops the attribute, so a subset is always regrouped.
  subset <- states[c(1L, 2L)]
  expect_null(attr(subset, "iwmde_active_groups", exact = TRUE))
  expect_identical(attr(.iwmde_row_state_groups(context, subset),
                        "row_edges", exact = TRUE), c(3L, 7L))

  # An empty list has no edges and still groups.
  empty <- .iwmde_row_state_groups(context, list())
  expect_identical(attr(empty, "n_states", exact = TRUE), 0L)
  expect_null(attr(empty, "row_edges", exact = TRUE))
})
