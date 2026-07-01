#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 14:56:49 2025

@author: valentin
"""
#!/home/valentin/anaconda3/envs/myenv python3
import cv2
import numpy as np
import os
import shutil

H_RANGE = (47, 49)
S_RANGE = (219, 221)
V_RANGE = (230, 232)

def calculate_color_percentage(image_path):
    image = cv2.imread(image_path)
    if image is None:
        print(f"[!] Could not read image: {image_path}")
        return 0.0  # Ensure it always returns a float
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    
    mask = (
        (hsv[:, :, 0] >= H_RANGE[0]) & (hsv[:, :, 0] <= H_RANGE[1]) &
        (hsv[:, :, 1] >= S_RANGE[0]) & (hsv[:, :, 1] <= S_RANGE[1]) &
        (hsv[:, :, 2] >= V_RANGE[0]) & (hsv[:, :, 2] <= V_RANGE[1])
    )
    return (np.count_nonzero(mask) / mask.size) * 100

def list_images_with_high_color_percentage(folder_path, threshold):
    high_color_images = []
    for filename in os.listdir(folder_path):
        if not filename.lower().endswith(".png"):
            continue
        image_path = os.path.join(folder_path, filename)
        color_pct = calculate_color_percentage(image_path)
        if color_pct > threshold:
            high_color_images.append(filename)
            print(f"[✓] {filename} - Color: {color_pct:.2f}%")
    return high_color_images

def move_high_color_images(src_folder, dst_folder, image_list):
    os.makedirs(dst_folder, exist_ok=True)
    for filename in image_list:
        src = os.path.join(src_folder, filename)
        dst = os.path.join(dst_folder, filename)
        if os.path.exists(src):
            shutil.move(src, dst)
            print(f"Moved: {filename}")
        else:
            print(f"[!] Missing: {src}")

def process_folders(base_dir, threshold):
    for subfolder in os.listdir(base_dir):
        path = os.path.join(base_dir, subfolder)
        #if not os.path.isdir(path): continue
        labelled = os.path.join(path, "labelled_tiles")
        original = os.path.join(path, "original_tiles")
        selected = os.path.join(path, "selected_tiles")
        if os.path.isdir(labelled) and os.path.isdir(original):
            print(f"\nProcessing: {subfolder}")
            high_color_imgs = list_images_with_high_color_percentage(labelled, threshold)
            move_high_color_images(original, selected, high_color_imgs)
        else:
            print(f"[!] Skipped: {subfolder} (missing folders)")

if __name__ == "__main__":
    base_path = "/home/valentin/Desktop/hepatovac/tiles/"
    threshold = 90
    process_folders(base_path, threshold)
