#!/bin/bash

# 脚本用途：循环执行 trento 命令生成数据，然后将数据转换为 ROOT 文件。
# 用法：./run_and_convert.sh <input_hdf5_path> <output_dir_path> <loop_count>
# 例如：./run_and_convert.sh /path/to/wosrc.hdf5 /path/to/output/dir 10

# 检查参数数量
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_hdf5_path> <output_dir_path> <loop_count>"
    echo "Example: $0 /Users/lipei/data/wosrc.hdf5 /Users/lipei/results 10"
    exit 1
fi

# 从命令行参数获取值
INPUT_HDF5_PATH=$1
OUTPUT_DIR=$2
N=$3

# 从输入文件路径中提取文件名作为基础名称
INPUT_FILENAME=$(basename "$INPUT_HDF5_PATH")
INPUT_BASE="${INPUT_FILENAME%.*}"

# 确保输入文件和输出目录存在
if [ ! -f "${INPUT_HDF5_PATH}" ]; then
    echo "Error: Input file ${INPUT_HDF5_PATH} not found!"
    exit 1
fi

if [ ! -d "${OUTPUT_DIR}" ]; then
    echo "Error: Output directory ${OUTPUT_DIR} not found!"
    exit 1
fi

# 定义基础 trento 命令
TRENTO_CMD="trento ${INPUT_HDF5_PATH} ${INPUT_HDF5_PATH} 1000000 -n 20 -x 7.09 -w 0.4 -m 6 -p 0.0 -k 100000 -v 0.4"

# ==========================================
# 阶段 1: 生成 .dat 文件
# ==========================================
echo "--- Starting Stage 1: Generating .dat files ---"
for ((i=1; i<=N; i++))
do
    OUTPUT_DAT_FILE="${OUTPUT_DIR}/${INPUT_BASE}_${i}.dat"
    echo "Running trento for: ${OUTPUT_DAT_FILE}"
    ${TRENTO_CMD} > "${OUTPUT_DAT_FILE}" 2>&1
    if [ $? -eq 0 ]; then
        echo "Successfully created ${OUTPUT_DAT_FILE}"
    else
        echo "Error: trento command failed for ${OUTPUT_DAT_FILE}"
    fi
done

# ==========================================
# 阶段 2: 将 .dat 文件转换为 .root 文件
# ==========================================
echo "--- Starting Stage 2: Converting .dat files to .root files ---"
for ((i=1; i<=N; i++))
do
    INPUT_DAT_FILE="${OUTPUT_DIR}/${INPUT_BASE}_${i}.dat"
    OUTPUT_ROOT_FILE="${OUTPUT_DIR}/${INPUT_BASE}_${i}.root"

    if [ ! -f "${INPUT_DAT_FILE}" ]; then
        echo "Warning: .dat file not found for conversion: ${INPUT_DAT_FILE}"
        continue
    fi

    echo "Converting ${INPUT_DAT_FILE} to ${OUTPUT_ROOT_FILE}"

    ./TransToRoot "${INPUT_DAT_FILE}" "${OUTPUT_ROOT_FILE}"

    if [ $? -eq 0 ]; then
        echo "Successfully converted to ${OUTPUT_ROOT_FILE}"
    else
        echo "Error: Conversion failed for ${INPUT_DAT_FILE}"
    fi
done

# ==========================================
# 阶段 3: 合并 .root 文件
# ==========================================
echo "--- 阶段 3: 正在合并 .root 文件 ---"
if ls "${OUTPUT_DIR}/${INPUT_BASE}_"*.root 1> /dev/null 2>&1; then
    MERGED_ROOT_FILE="${OUTPUT_DIR}/${INPUT_BASE}.root"
    echo "正在将所有 .root 文件合并到 ${MERGED_ROOT_FILE}"
    hadd -f "${MERGED_ROOT_FILE}" "${OUTPUT_DIR}/${INPUT_BASE}_"*.root
    if [ $? -eq 0 ]; then
        echo "合并成功！合并后的文件: ${MERGED_ROOT_FILE}"

        # ==========================================
        # 阶段 4: 清理临时文件
        # ==========================================
        echo "--- 阶段 4: 正在清理临时文件 ---"
        echo "正在删除所有临时 .dat 和 .root 文件..."
        rm "${OUTPUT_DIR}/${INPUT_BASE}_"*.dat "${OUTPUT_DIR}/${INPUT_BASE}_"*.root
        echo "清理完成。"
    else
        echo "错误: hadd 合并失败。未进行清理。"
    fi
else
    echo "警告: 没有找到 .root 文件进行合并。"
fi


echo "Script finished."