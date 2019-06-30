# GTF-Compare
## Description
This tool is used to compare GTF files with a reference GTF.
Currently it generates 6 files for each query GTF as follows:
<table>
    <thead>
        <th>File Name</th>
        <th>File Extension</th>
        <th>Description</th>
    </thead>
    <tbody>
        <tr>
            <td rowspan="3">&lt;query GTF&gt;</td>
            <td>exon</td>
            <td>Reports the best match (if found) for each exon found in the query GTF</td>
        </tr>
        <tr>
            <td>itracking</td>
            <td>Reports the query GTF interval structure and it's overlap with reference GTF if found</td>
        </tr>
        <tr>
            <td>tstat</td>
            <td>Temporary simple statistics file on the query GTF to be used later within the tool</td>
        </tr>
        <tr>
            <td rowspan="2">&lt;ref_query GTF&gt;</td>
            <td>exon</td>
            <td>Reports the best match (if found) for each exon found in the reference GTF</td>
        </tr>
        <tr>
            <td>itracking</td>
            <td>Reports the reference GTF interval structure and it's overlap with the query GTF if found</td>
        </tr>
        <tr>
            <td>ref</td>
            <td>tstat</td>
            <td>Temporary simple statistics file on the reference GTF to be used later within the tool</td>
        </tr>
    </tbody>
</table>

---
## How to install
### Install prerequisites
```
# Install pypy
sudo apt install pypy

# Download pip
wget https://bootstrap.pypa.io/get-pip.py
pypy get-pip.py

# Install intervaltree
pypy -m pip install intervaltree
pypy -m pip install numpy
```
### Get the source code
```
git clone https://github.com/abdelrahmanMA/gtf-compare.git
```
---
### How to use

```
cd gtf-compare/gtfcompare
pypy compare.py -r [<reference.gtf>] [<query.gtf> [<query.gtf> ...]]
```
You can pass as many query GTFs as you want.
Alternatively -i can be used to pass a text file with a list of query GTFs

```
pypy compare.py -r [<reference.gtf>] -i [<list.txt>]
```
This tool can be used with python 2.7 and python 3.6, but for the sake of speed pypy is recommended.