# 导入所需的库
import numpy as np
import h5py
import os
import glob

def convert_matlab_array_to_hdf5(input_filename, output_filename):
    """
    此函数将从MATLAB导出的数据文件转换为HDF5文件。
    
    参数:
    input_filename (str): 从MATLAB导出的数据文件路径（例如 'nuclei_coords.txt'）。
    output_filename (str): 输出的HDF5文件路径（例如 'nuclei_data.hdf5'）。
    """
    
    print(f"正在读取文件: {input_filename}...")
    
    # 检查输入文件是否存在
    if not os.path.exists(input_filename):
        print(f"错误: 输入文件 '{input_filename}' 未找到。请确保它存在于当前目录下。")
        return
        
    try:
        # 使用NumPy读取数据。dtype=np.float64是为了确保精度与MATLAB的双精度一致。
        data_2d = np.loadtxt(input_filename, dtype=np.float64)
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return

    # 打印原始数组的形状以进行验证
    print(f"读取的原始数组形状: {data_2d.shape}")
    
    # 获取Ns（核的数量）
    Ns = data_2d.shape[0]
    
    # 获取每个核子的坐标总数，即48
    num_coords_per_nucleus = data_2d.shape[1]
    
    # 检查维度是否正确（Ns x 48）
    if num_coords_per_nucleus != 36:
        print("错误: 输入数组的第二维不是48。请检查你的MATLAB数据。")
        return
        
    # 重塑数组到所需的形状 (Ns, 48, 3)
    # 这将把48个坐标（x1,y1,z1,...x16,y16,z16）重新组织成16个核子，每个核子有3个坐标。
    data_3d = data_2d.reshape(Ns, 12, 3)
    
    # 打印重塑后的数组形状以进行验证
    print(f"重塑后的数组形状: {data_3d.shape}")
    
    # 使用h5py创建HDF5文件
    try:
        with h5py.File(output_filename, 'w') as f:
            # 在HDF5文件中创建一个名为 'nuclei_coordinates' 的数据集
            # 并将重塑后的数据写入其中
            dset = f.create_dataset('nuclei_coordinates', data=data_3d)
            print(f"数据已成功保存到HDF5文件: {output_filename}")
            print(f"HDF5文件中的数据集名称: 'nuclei_coordinates'")
            print(f"HDF5数据集的形状: {dset.shape}")
            
    except Exception as e:
        print(f"写入HDF5文件时发生错误: {e}")

# --- 如何使用此脚本 ---
# 自动查找当前目录下的唯一 .txt 文件，然后进行转换
if __name__ == "__main__":
    # 获取当前目录下所有 .txt 文件
    txt_files = glob.glob("*.txt")
    
    # 根据实际情况选择文件
    if len(txt_files) == 0:
        print("错误: 当前目录下没有找到任何 .txt 文件。")
        exit(1)
    elif len(txt_files) > 1:
        print("错误: 当前目录下存在多个 .txt 文件，请确保只有一个。找到的文件列表如下:")
        for f in txt_files:
            print(f"  - {f}")
        exit(1)
    else:
        input_file = txt_files[0]
        print(f"自动找到输入文件: {input_file}")
    
    output_file = 'initial.hdf5'
    convert_matlab_array_to_hdf5(input_file, output_file)

