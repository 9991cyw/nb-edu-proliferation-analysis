# ============================================================
# 批量在原始图像上画红色虚线框（整图）并标注编号
# 用于投稿时提交"原始数据（not for review）"材料
#
# 用法：
#   1) 准备清单 CSV（两列：file, label），UTF-8 编码，无表头：
#      IMG_0231.tif,Fig 4B panel1
#      IMG_0237.tif,Fig 4B panel2
#   2) 修改下方 input_dir / output_dir / manifest 路径
#   3) 运行：Rscript draw_dashed_box.R
# ============================================================
library(magick)

# ---- 配置区（按需修改）----
# 建议自己新建三个位置，脚本只读取，绝不改动原始文件：
#   D:/标注任务/原始图    -> 放被选中的原始图
#   D:/标注任务/清单.csv  -> 标注清单
#   D:/标注任务/标注版    -> 输出（可不存在，脚本自动创建）
input_dir  <- "D:/标注任务/原始图"   # 1. 被选中原始图所在文件夹
output_dir <- "D:/标注任务/标注版"   # 2. 标注版输出文件夹
manifest   <- "D:/标注任务/清单.csv" # 3. 清单路径；设为 "" 则自动处理 input_dir 内全部图片
# 清单 CSV 格式（用 Excel 另存为"CSV UTF-8"）：
#   第一列 = 图片文件名，必须与 input_dir 里的文件名完全一致（含扩展名、大小写）
#   第二列 = 显示在图片右上角的文字，例如：
#     NB1+SC_0014.tif,Figure 5A panel1
#     NB1+SC_0021.tif,Figure 5A panel2
#     NB2+SC_0004.tif,Figure 5B panel1
box_color  <- "red"                       # 虚线框颜色
box_lwd    <- 6                           # 虚线框线宽（像素）
label_cex  <- 3                           # 编号字号
label_color<- "red"                       # 编号颜色

# ---- 以下一般无需修改 ----
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 图片格式白名单
ext_ok <- tolower(c("tif","tiff","png","jpg","jpeg","bmp","gif"))

if (nzchar(manifest) && file.exists(manifest)) {
  m <- read.csv(manifest, header = FALSE, stringsAsFactors = FALSE,
                col.names = c("file", "label"),
                fileEncoding = "UTF-8-BOM")   # Excel 保存的 CSV UTF-8 带 BOM，需特殊处理
} else {
  m <- data.frame(
    file = list.files(input_dir, full.names = FALSE),
    label = list.files(input_dir, full.names = FALSE),
    stringsAsFactors = FALSE
  )
}
m <- m[tolower(tools::file_ext(m$file)) %in% ext_ok, ]
if (nrow(m) == 0) stop("没有找到可处理的图片，请检查路径与格式")

for (i in seq_len(nrow(m))) {
  src <- file.path(input_dir, m$file[i])
  img <- tryCatch(image_read(src), error = function(e) {
    cat("跳过(无法读取):", m$file[i], "\n"); NULL
  })
  if (is.null(img)) next
  info <- image_info(img)
  W <- info$width; H <- info$height

  # 用 R 图形设备叠加：点线框 + 编号（期刊原文要求 dotted line 标记）
  img2 <- image_draw(img)
  rect(1, 1, W - 1, H - 1,
       border = box_color, lwd = box_lwd, lty = 3)          # lty=3 点线，框住整图
  text(W - 10, H - 10, labels = m$label[i],
       adj = c(1, 1), col = label_color, cex = label_cex)    # 右上角编号
  dev.off()

  stem  <- tools::file_path_sans_ext(m$file[i])
  ext   <- tools::file_ext(m$file[i])
  out   <- file.path(output_dir, paste0(stem, "_marked.", ext))
  image_write(img2, out)
  cat("已处理:", m$file[i], "->", basename(out), " (", W, "x", H, ")\n")
}
cat("全部完成，共", nrow(m), "张\n")
