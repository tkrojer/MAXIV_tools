import sqlalchemy
from sqlalchemy.sql import select
from sqlalchemy import and_, or_

import sys
sys.path.append('/data/staff/biomax/tobias/software/MAXIV_tools/lib')
from db import dal

import pandas as pd
from fpdf import FPDF
import matplotlib
matplotlib.use('Agg')  # Use the Agg backend, otherwise the backend does rely on qt which might not be set up correctly
from matplotlib import image as mpimg
from matplotlib import pyplot as plt

db_file = "/data/visitors/biomax/20240919/20240317/fragmax/lab/database/fragmax.sqlite"
tmp_dir = "/data/visitors/biomax/20240919/20240317/fragmax/tmp"
dal.db_init(db_file)

k = dal.soaked_crystals_table.join(
        dal.marked_crystals_table, dal.soaked_crystals_table.c.marked_crystal_id ==
                                   dal.marked_crystals_table.c.marked_crystal_id, isouter=True).join(
        dal.mounted_crystals_table, dal.marked_crystals_table.c.marked_crystal_id ==
                                    dal.mounted_crystals_table.c.marked_crystal_id, isouter=True).join(
        dal.xray_dataset_table, dal.mounted_crystals_table.c.mounted_crystal_id ==
                                dal.xray_dataset_table.c.mounted_crystal_id, isouter=True).join(
        dal.xray_processing_table, dal.xray_dataset_table.c.dataset_id ==
                                   dal.xray_processing_table.c.dataset_id, isouter=True).join(
        dal.xray_initial_refinement_table, dal.xray_processing_table.c.processing_id ==
                                           dal.xray_initial_refinement_table.c.processing_id, isouter=True).join(
        dal.soak_plate_table, dal.soaked_crystals_table.c.soak_plate_id ==
                              dal.soak_plate_table.c.soak_plate_id, isouter=True).join(
        dal.compound_batch_table, dal.soak_plate_table.c.compound_batch_code ==
                                  dal.compound_batch_table.c.compound_batch_code, isouter=True).join(
        dal.compound_table, dal.compound_batch_table.c.compound_code ==
                            dal.compound_table.c.compound_code, isouter=True).join(
        dal.crystal_screen_condition_table, dal.marked_crystals_table.c.crystal_screen_condition_id ==
                                            dal.crystal_screen_condition_table.c.crystal_screen_condition_id,
        isouter=True).join(
        dal.crystal_plate_table, dal.marked_crystals_table.c.crystal_plate_id ==
                                 dal.crystal_plate_table.c.crystal_plate_id, isouter=True).join(
        dal.protein_batch_table, dal.crystal_plate_table.c.protein_batch_id ==
                                 dal.protein_batch_table.c.protein_batch_id, isouter=True)

q = select([dal.mounted_crystals_table.c.mounted_crystal_code,
            dal.soak_plate_table.c.base_buffer,
            dal.marked_crystals_table.c.marked_crystal_image,
            dal.xray_dataset_table.c.crystal_snapshot_1,
            dal.xray_dataset_table.c.crystal_snapshot_2,
            dal.xray_dataset_table.c.crystal_snapshot_3,
            dal.xray_dataset_table.c.crystal_snapshot_4,
            dal.xray_dataset_table.c.is_dataset,
            dal.xray_dataset_table.c.selected.label('xray_dataset_table_selected'),
            dal.xray_processing_table.c.selected.label('xray_processing_table_selected'),
            dal.xray_processing_table.c.reflns_d_resolution_high,
            dal.xray_processing_table.c.reflns_outer_pdbx_netI_over_sigmaI,
            dal.xray_processing_table.c.reflns_inner_pdbx_Rmerge_I_obs,
            dal.xray_processing_table.c.sym_space_group
           ]).order_by(dal.mounted_crystals_table.c.mounted_crystal_code.asc())
q = q.select_from(k)
df = pd.read_sql_query(q, dal.engine)
df = df.dropna()

df = df[df["is_dataset"] != False]
df = df[df["xray_dataset_table_selected"] != False]
df = df[df["xray_processing_table_selected"] != False]

df = df.drop('is_dataset', axis=1)
df = df.drop('xray_dataset_table_selected', axis=1)
df = df.drop('xray_processing_table_selected', axis=1)

df.rename(columns={'mounted_crystal_code': 'sample_id',
                   'crystal_snapshot_1': 'img1',
                   'crystal_snapshot_2': 'img2',
                   'crystal_snapshot_3': 'img3',
                   'crystal_snapshot_4': 'img4',
                   'reflns_d_resolution_high': 'reso_high',
                   'reflns_outer_pdbx_netI_over_sigmaI': 'I_sigI',
                   'sym_space_group': 'spg',}, inplace=True)


class PDF(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 12)

    #        self.cell(0, 10, 'Images and Text in Same Row', 0, 1, 'C')

    def footer(self):
        #        self.set_y(-15)
        self.set_font('Arial', 'I', 8)


#        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

# pdf = PDF()
pdf = PDF(format=(1000, 1400))  # Custom page size: width=210mm, height=400mm
pdf.add_page()
pdf.set_font('Arial', '', 12)

image_width_mm = 60  # Image width in mm
image_height_mm = 40  # Image height in mm
column_width = 60  # Column width for text in mm

# Table header
column_widths = [30, 140, 60, 60, 60, 60, 60, 30, 30, 30]  # Adjust the widths as necessary
pdf.set_font('Arial', 'B', 12)
headers = ["sample_id", "base_buffer", "marked_crystal_image", "img1", "img2", "img3", "img4", "reso", "i/sig(I)", "spg"]
for i, header in enumerate(headers):
    pdf.cell(column_widths[i], 10, header, 1, 0, 'C')
pdf.ln(10)

for index, row in df.iterrows():
    print(f">>>> index: {index}")
    # Determine start position of the image
    #    image_x = pdf.get_x() + column_width * 2 + 100
    image_x = pdf.get_x() + column_widths[0] + column_widths[1]
    image_a = pdf.get_x() + column_widths[0] + column_widths[1] + 1 * image_width_mm
    image_b = pdf.get_x() + column_widths[0] + column_widths[1] + 2 * image_width_mm
    image_c = pdf.get_x() + column_widths[0] + column_widths[1] + 3 * image_width_mm
    image_d = pdf.get_x() + column_widths[0] + column_widths[1] + 4 * image_width_mm
    image_y = pdf.get_y()

    # Add text cells
    pdf.cell(column_widths[0], image_height_mm, str(row['sample_id']), border=1)
    pdf.cell(column_widths[1], image_height_mm, str(row['base_buffer']), border=1)

    # Add the image
    img_path = row['marked_crystal_image'].replace('/Users/tobkro/tmp/20240111', '/data/visitors/biomax/20240919/20240317/fragmax/lab')
    if img_path:  # Check if the image path is not empty or None
        img = mpimg.imread(img_path)
        plt.imshow(img)
        plt.axis('off')  # Do not display axis
        temp_image_path = f'{tmp_dir}/temp_image_{index}.png'
        print(img_path, temp_image_path)
        plt.savefig(temp_image_path, bbox_inches='tight', pad_inches=0)
        plt.close()
        pdf.image(temp_image_path, x=image_x, y=image_y, w=image_width_mm, h=image_height_mm)

    img_path = row['img1']
    if img_path:  # Check if the image path is not empty or None
        img = mpimg.imread(img_path)
        plt.imshow(img)
        plt.axis('off')  # Do not display axis
        temp_image_path = f'{tmp_dir}/temp_image_{index}_1.png'
        print(img_path, temp_image_path)
        plt.savefig(temp_image_path, bbox_inches='tight', pad_inches=0)
        plt.close()
        pdf.image(temp_image_path, x=image_a, y=image_y, w=image_width_mm, h=image_height_mm)

    img_path = row['img2']
    if img_path:  # Check if the image path is not empty or None
        img = mpimg.imread(img_path)
        plt.imshow(img)
        plt.axis('off')  # Do not display axis
        temp_image_path = f'{tmp_dir}/temp_image_{index}_2.png'
        print(img_path, temp_image_path)
        plt.savefig(temp_image_path, bbox_inches='tight', pad_inches=0)
        plt.close()
        pdf.image(temp_image_path, x=image_b, y=image_y, w=image_width_mm, h=image_height_mm)

    img_path = row['img3']
    if img_path:  # Check if the image path is not empty or None
        img = mpimg.imread(img_path)
        plt.imshow(img)
        plt.axis('off')  # Do not display axis
        temp_image_path = f'{tmp_dir}/temp_image_{index}_3.png'
        print(img_path, temp_image_path)
        plt.savefig(temp_image_path, bbox_inches='tight', pad_inches=0)
        plt.close()
        pdf.image(temp_image_path, x=image_c, y=image_y, w=image_width_mm, h=image_height_mm)

    img_path = row['img4']
    if img_path:  # Check if the image path is not empty or None
        img = mpimg.imread(img_path)
        plt.imshow(img)
        plt.axis('off')  # Do not display axis
        temp_image_path = f'{tmp_dir}/temp_image_{index}_4.png'
        print(img_path, temp_image_path)
        plt.savefig(temp_image_path, bbox_inches='tight', pad_inches=0)
        plt.close()
        pdf.image(temp_image_path, x=image_d, y=image_y, w=image_width_mm, h=image_height_mm)

    # Add text cells
    pdf.cell(column_widths[7], image_height_mm, str(row['reso_high']), border=1)
    pdf.cell(column_widths[8], image_height_mm, str(row['I_sigI']), border=1)
    pdf.cell(column_widths[9], image_height_mm, str(row['spg']), border=1)
#    pdf.cell(column_widths[7], image_height_mm, str(row['reso_high']), border=1)
#    pdf.cell(column_widths[8], image_height_mm, str(row['I_sigI']), border=1)
#    pdf.cell(column_widths[9], image_height_mm, str(row['spg']), border=1)

    # Move below the image/row for the next entry
    pdf.ln(image_height_mm)

pdf.output('images_text_same_row.pdf')
