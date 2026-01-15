# scripts/01_hello_boundary_logic.R
# Purpose: sanity check that the repo runs and produces an output.
# test merge
message("Hello from boundary_logic (feature branch).")

set.seed(1)
df <- data.frame(
  x = rnorm(200),
  y = rnorm(200)
)
df$class <- ifelse(df$x + df$y > 0, "A", "B")

# Create output directory if needed
if (!dir.exists("outputs")) dir.create("outputs")

# Plot with a decision boundary
png(filename = "outputs/hello_scatter.png", width = 900, height = 700)
plot(df$x, df$y,
     col = ifelse(df$class == "A", "blue", "red"),
     pch = 19,
     main = "boundary_logic: first decision boundary",
     xlab = "x",
     ylab = "y")

abline(a = 0, b = -1, lwd = 2)  # x + y = 0
legend("topright", legend = c("Class A", "Class B"),
       col = c("blue", "red"), pch = 19)
dev.off()

message("Wrote outputs/hello_scatter.png with decision boundary")
