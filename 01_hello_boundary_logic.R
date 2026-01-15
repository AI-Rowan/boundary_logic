# scripts/01_hello_boundary_logic.R
# Purpose: sanity check that the repo runs and produces an output.

message("Hello from boundary_logic.")

# Minimal example dataset
set.seed(1)
df <- data.frame(
  x = rnorm(200),
  y = rnorm(200)
)
df$class <- ifelse(df$x + df$y > 0, "A", "B")

# Simple plot output (writes a file so you see reproducibility)
if (!dir.exists("outputs")) dir.create("outputs")

png(filename = "outputs/hello_scatter.png", width = 900, height = 700)
plot(df$x, df$y, main = "boundary_logic: first output",
     xlab = "x", ylab = "y")
dev.off()

message("Wrote outputs/hello_scatter.png")
