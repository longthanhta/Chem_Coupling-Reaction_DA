import fitz # PyMuPDF
import io
from PIL import Image
import os

file = "ja044588b.pdf"
try:
    os.mkdir('figures')
except:
    print('figures folder already exis')
parent_path=os.getcwd()
for path, subdirs, files in os.walk(os.getcwd()):
    for filename in files:
        if '.pdf' not in filename: continue
        pdf_file = fitz.open(filename)
        try:
            filename=filename.replace('.pdf','')
            fig_path='figures/'+filename
            os.mkdir(fig_path)
        except:
            print(fig_path+'folder already exis')
        # iterate over PDF pages
        for page_index in range(len(pdf_file)):
            # get the page itself
            page = pdf_file[page_index]
            image_list = page.getImageList()
            # printing number of images found in this page
            if image_list:
                print(f"[+] Found a total of {len(image_list)} images in page {page_index}")
            else:
                print("[!] No images found on page", page_index)
            for image_index, img in enumerate(page.getImageList(), start=1):
                # get the XREF of the image
                xref = img[0]
                # extract the image bytes
                base_image = pdf_file.extractImage(xref)
                image_bytes = base_image["image"]
                # get the image extension
                image_ext = base_image["ext"]
                # load it to PIL
                image = Image.open(io.BytesIO(image_bytes))
                # save it to local disk
                image_filename=f"image{page_index+1}_{image_index}.{image_ext}"
                print('fig_path',fig_path)
                image.save(open(fig_path+'/'+image_filename, "wb"))
