# 遍历所有以 _1.clean.fq.gz 或 _2.clean.fq.gz 结尾的文件
for file in *_[12].clean.fq.gz; do
    # 提取 SampleID（第一个下划线之前的部分）
    sample_id=$(echo "$file" | cut -d'_' -f1)
    # 提取文件后缀（_1.clean.fq.gz 或 _2.clean.fq.gz）
    suffix=$(echo "$file" | grep -o '_[12].clean.fq.gz')

    # 组合新的文件名
    new_name="${sample_id}${suffix}"

    # 输出原文件名和新文件名
    echo "Renaming: $file -> $new_name"

    # 实际进行重命名
    mv "$file" "$new_name"
done
