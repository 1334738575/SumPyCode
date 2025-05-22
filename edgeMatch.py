import cv2
import numpy as np

def showImg(img, winName, waitKey=False):
    cv2.imshow(winName, img)
    if waitKey:
        cv2.waitKey(0)

debugShow = True

# 读取图像和模板
raw_image = cv2.imread('C:/Users/DELL/Desktop/code/GeoMatch_demo/Search1.jpg', cv2.IMREAD_GRAYSCALE)
raw_template = cv2.imread('C:/Users/DELL/Desktop/code/GeoMatch_demo/Template.jpg', cv2.IMREAD_GRAYSCALE)

# # 指定新的尺寸
# new_width = 800
# new_height = 600
# new_size = (new_width, new_height)
# # 使用 cv2.resize() 进行缩放
# image = cv2.resize(raw_image, new_size)
# template = cv2.resize(raw_template, new_size)
image = raw_image
template = raw_template

# # 显示原始图像和缩放后的图像
# cv2.imshow('Original Image', image)
# cv2.imshow('Resized Image', resized_image)
# # 等待键盘输入并关闭窗口
# cv2.waitKey(0)
# cv2.destroyAllWindows()
# exit(0)


# 边缘检测
image_edges = cv2.Canny(image, threshold1=50, threshold2=150)
template_edges = cv2.Canny(template, threshold1=50, threshold2=150)
if debugShow:
    showImg(image_edges, "edge1", False)
    showImg(template_edges, "edge2", True)

# 模板匹配
result = cv2.matchTemplate(image_edges, template_edges, cv2.TM_CCOEFF_NORMED)

# 找到最佳匹配位置
min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(result)

# 获取模板的宽度和高度
h, w = template_edges.shape

# 绘制矩形框标记匹配位置
top_left = max_loc
bottom_right = (top_left[0] + w, top_left[1] + h)
cv2.rectangle(image, top_left, bottom_right, 255, 2)

# 显示结果
cv2.imshow('Matched Result', image)
cv2.waitKey(0)
cv2.destroyAllWindows()
