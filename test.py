from spmUtils import *

#base_path =f"/media/captainbroccoli/DATA"
base_path = f"/media/captainbroccoli/DATA"
#folder="2022-07-17"
folder="2022-07-27"
image_prename = f"20220727-154751_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--"

for i in range(342):
    image_name= f"{image_prename}{i+1}_1.Z_mtrx"
    #image_name= "20220729-095353_Cu(111)--AFM_NonContact_QPlus_AtomManipulation_AuxChannels--168_1.Z_mtrx"
    image_path = f"{base_path}/{folder}/{image_name}"
    image, message, message_im= get_matrix_image(image_path)

    # plt.imshow(image.data, interpolation='nearest')
    # plt.show()




