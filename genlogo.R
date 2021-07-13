#install.packages("hexSticker")
library(hexSticker)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("emojifont")
library(emojifont)
# fa4
#r <- ggplot() + geom_fontawesome("fas fa-shoe-prints", color='yellow', size = 20) + theme_void()

# fa5
#install.packages("fontawesome")
library(fontawesome)
aaa=fa(name = "fas fa-shoe-prints", fill = "#F6F622")

print_fontawesome<- function (x, view = interactive(), ...)
{
  dots <- list(...)
  html <- paste(x, collapse = "\n")
  c("<!DOCTYPE html>", "<html>", "<head>", "<meta charset=\"utf-8\">",
    "</head>", "<body>", html, "</body>", "</html>") %>%
    paste(collapse = "\n") %>% htmltools::HTML() %>% htmltools::html_print()
  return(html)
}

library(rsvg)
cat(print_fontawesome(aaa, view=FALSE),file="man/figures/fa-shoe-prints.svg")
SVGres <- rsvg::rsvg_svg("man/figures/fa-shoe-prints.svg","man/figures/Cairo-fa-shoe-prints.svg")
SVGres <- rsvg::rsvg_svg("man/figures/fa-shoe-prints.svg")
PDFres <- rsvg::rsvg_pdf("man/figures/fa-shoe-prints.svg","man/figures/Cairo-fa-shoe-prints.pdf")
PDFres <- rsvg::rsvg_pdf("man/figures/fa-shoe-prints.svg")

library(png)
library(grid)
library(cowplot)
library(magick)
theme_set(theme_cowplot())

r = ggdraw() +
  draw_image("man/figures/Cairo-fa-shoe-prints.pdf", scale = .55)

sticker(
  r,
  package = "bootPLS",
  p_size = 7,
  s_x = .98,
  s_y = 0.7,
  s_width = 1.7,
  s_height = 1.3,
  p_x = 1,
  p_y = 1.3,
  url = "https://cran.r-project.org/package=bootPLS",
  u_color = "white",
  u_size = 1.05,
  h_fill = "black",
  h_color = "grey",
  filename = "man/figures/logo.png"
)
