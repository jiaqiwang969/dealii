#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)#!/usr/bin/env python
#-----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
#-----------------------------------------
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.h$',args.filename)==None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.h$','',args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()



text=source

## 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r"^\s\*|^\/\/.*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)

output_filename = filebase+'_pre.h'
with open(output_filename, 'w') as txt_file:
    txt_file.write(text)
