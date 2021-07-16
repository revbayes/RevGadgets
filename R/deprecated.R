.add_epoch_times <- function(p, max_age, dy_bars, dy_text) {
  max_x <- max(p$data$x)
  max_y <- max(p$data$y)

  epoch_names <- rev(
    c(
      "Holocene",
      "Pleistocene",
      "Pliocene",
      "Miocene",
      "Oligocene",
      "Eocene",
      "Paleocene",
      "Upper\nCretaceous",
      "Lower\nCretaceous"
    )
  )
  epoch_ages <- rev(c(0, 0.117, 2.588, 5.332, 23.03,
                      33.9, 55.8, 65.5, 99.6))
  period_names <- rev(
    c(
      "Quaternary",
      "Neogene",
      "Paleogene",
      "Cretaceous",
      "Jurassic",
      "Triassic",
      "Permian",
      "Carboniferous",
      "Devonian",
      "Silurian",
      "Ordovician",
      "Cambrian"
    )
  )

  period_ages <- rev(c(
    0,
    2.588,
    23.03,
    65.5,
    145.5,
    199.6,
    251,
    299,
    359.2,
    416,
    443.7,
    488.3
  ))
  if (max_age > 140) {
    x_pos <- max_x - c(max_age, period_ages)
  } else
    x_pos <- max_x - c(max_age, epoch_ages)

  y_pos <- rep(max_y, length(x_pos))
  x_pos_mid <-
    (x_pos[1:(length(x_pos) - 1)] + x_pos[2:length(x_pos)]) / 2

  for (k in 2:(length(x_pos))) {
    box_col <- "gray92"
    if (k %% 2 == 0)
      box_col <- "white"
    box <-
      ggplot2::geom_rect(
        xmin = x_pos[k - 1],
        xmax = x_pos[k],
        ymin = dy_bars,
        ymax = y_pos[k],
        fill = box_col
      )
    p <- gginnards::append_layers(p, box, position = "bottom")
  }
  if (max_age > 140) {
    for (k in seq_len(length(period_names))) {
      p <-
        p + ggplot2::annotate(
          geom = "text",
          label = period_names[k],
          angle = 90,
          x = x_pos_mid[k],
          y = dy_text,
          hjust = 0,
          size = 3.25
        )
    }
  } else {
    for (k in seq_len(length(epoch_names))) {
      p <-
        p + ggplot2::annotate(
          geom = "text",
          label = epoch_names[k],
          angle = 90,
          x = x_pos_mid[k],
          y = dy_text,
          hjust = 0,
          size = 3.25
        )
    }
  }
  return(p)
}

# stolen from treeio: https://github.com/YuLab-SMU/treeio
.validTblTree <-
  function(object, cols = c("parent", "node", "label")) {
    cc <- cols[!cols %in% colnames(object)]
    if (length(cc) > 0) {
      msg <-
        paste0("invalid tbl_tree object.\n  missing column:\n    ",
               paste(cc, collapse = ","),
               ".")
    }
  }
