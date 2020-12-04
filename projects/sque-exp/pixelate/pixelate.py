from PIL import Image

# choose picture
img = input("picture name : ")

# Open Paddington
img = Image.open(img)

# Resize smoothly down to 16x16 pixels
img_= img.resize((32,32),resample=Image.BILINEAR)

# Scale back up using NEAREST to original size
img = img_.resize(img.size,Image.NEAREST)

# Save
img.save('img.png')
