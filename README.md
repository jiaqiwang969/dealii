# deal.II-translator
Using Doxygen format to translate packages into different languages-https://www.dealii.org/

- 首先，我们利用"./to.py filename" 生成剔除特殊字符的中转文件；
- 然后，利用deepl或者google翻译；
- 然后，利用"./from.py filename_0_T.txt";
- 然后，得到的替换文件回代到source中，再一次进行编译。注意编译需要开启mathjax，否则公式可能乱码"cmake -DDEAL_II_COMPONENT_DOCUMENTATION=ON -DDEAL_II_DOXYGEN_USE_MATHJAX=ON",然后再输入"make documentation"，若有改动，只需输入后者，无需重头开始编译。

可能用到的批量操作代码：
- 拆分一次性翻译的文档：
awk '/\[2.x.0\]/{n++;w=1} n&&w{print >"step-"n"_0_T.txt"}' T.txt
- 重新批量命名
for file in `ls *.h`;do mv $file `echo $file|sed 's/_0_T//g'`;done;
- 批量允许python代码
for filename in headers/*.h; do
        ./to.py "$filename" 
done
for filename in examples/*_T.txt; do
        ./from.py ${filename} 
done
-在当前目录下，将所有aaaModule都替换为bbbName
grep -rl 'aaaModule' ./  | xargs sed -i "" "s/aaaModule/bbbName/g"

-r 表示搜索子目录
-l 表示输出匹配的文件名


# 将h文件全部翻译，并进行CI-测试代码 
使用说明：
- step1: terminal 运行 `translator_debug.sh > log1.txt` 用来测试，如果有错误，需要进行bug调试
- step2: terminal 运行 `translator_to.sh > log2.txt` 进行元素提取，整理打包文件在transltor_file 找到。 运行 `cat * > merge.txt` 合并为一个。
- step3: 在Testlab文件夹中，替换已经翻译好的merge.txt， 然后运行  `bash_T.sh` .
- step4: 最后进行整理，替换。 其中，tutorial有一些文件在修复wrapcomments过程中，出错为空，暂且替换为未修复版本。

# bug
- bug01: base/polunomials_p.h,fe/fe_simple_p.h 文件的代码块曾出现重复现象，，暂未解决，通过报错找到，并人工修复。 具体位置见commit。2021-06-07