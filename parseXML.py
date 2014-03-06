import xml.etree.cElementTree as ET

tree = ET.ElementTree(file='config/test_fullsamples.xml')
root = tree.getroot()

#root.tag, root.attrib

for child_of_root in root:
      child_of_root.set('add', '0')

#import sys
#tree.write(sys.stdout)


for rank in root.iter('add'):
#       new_rank = int(rank.text) + 1
#       rank.text = str(1)
#       rank.set('updated', 'yes')

        print rank.value
tree.write('output.xml')
