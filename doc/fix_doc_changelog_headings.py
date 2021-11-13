import os.path
import sys
"""This script corrects a problem in the file md_CHANGELOG.html, produced by make doc.
   Due to a bug in Doxygen, the links in a header line in a Markdown file are not
   properly converted in html. The symbols '<' and '>' that delimit the html tags
   '<a' and '/a>' are translated into the html codes '&lt;' and '&gt;', so a
   browser displays the html code of a link instead of a formatted link.
   This script simply replaces '&lt;' and '&gt;' with '<' and '>' if a line
   starts with an '&'. It has to be run after the compilation of the documentation.
   Once the bug in Doxygen is fixed, this script can be removed.
   See also the discussion in: https://github.com/smash-transport/smash-devel/issues/640.
"""

# We check that the script is invoked with the correct number of arguments,
# if not we print a reminder about the syntax
if (len(sys.argv)!=2):
    print("Syntax: fix_doc_changelog_headings.py <path/md_CHANGELOG.html>\n")
    sys.exit(1)

# We check that the file exists and that its name is what we expect to be
if ((not os.path.exists(sys.argv[1])) or (os.path.basename(sys.argv[1]) != "md_CHANGELOG.html")):
    print("Please, check that the argument of fixdoc.py corresponds to the path of the file md_CHANGELOG.html, \
            including the name of the file itself\n")
    sys.exit(1)

# We open the file and read all the lines, the result goes in the list 'lines'
with open(sys.argv[1],"r") as infile:
    lines=infile.readlines()

# In the lines that begin with the '&' character, we replace the html codes of '<' and '>'
# with these symbols
for i, line in enumerate(lines):
    if (line[0]=="&"):
        line_tmp=line.replace("&lt;","<")
        lines[i]=line_tmp.replace("&gt;",">")
        
# We overwrite the file md_CHANGELOG.html with the corrected lines
with open(sys.argv[1],"w") as outfile:
    outfile.writelines(lines)
