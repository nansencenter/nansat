#------------------------------------------------------------------------------
# Name:         test_node.py
# Purpose:      Test the Node class
#
# Author:       Aleksander Vines
#
# Created:      2016-02-26
# Last modified:2016-02-26T16:00
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------
from __future__ import absolute_import
import unittest
import os
from . import nansat_test_data as ntd
from nansat.node import Node


class NodeTest(unittest.TestCase):
    def test_creation(self):
        tag = 'Root'
        value = '   Value   '
        anAttr = 'elValue'
        new_value = 'New Value'
        node = Node(tag, value=value, anAttr=anAttr)
        self.assertEqual(node.tag, tag)
        self.assertDictEqual(node.attributes, {'anAttr': anAttr})
        self.assertEqual(node.value, value.strip())
        self.assertEqual(node[tag], value.strip())
        node[tag] = new_value
        self.assertEqual(node.value, new_value)

    def test_getAttributeList(self):
        tag = 'Root'
        value = '   Value   '
        anAttr = 'elValue'
        secondAttr = 'Some value'
        finalAttribute = 'A last value'
        node = Node(tag, value=value, anAttr=anAttr, secondAttr=secondAttr,
                    finalAttribute=finalAttribute)
        nameList, valList = node.getAttributeList()
        self.assertIsInstance(nameList, list)
        self.assertIsInstance(valList, list)
        index = valList.index(anAttr)
        self.assertEqual(nameList[index], 'anAttr')
        index = valList.index(secondAttr)
        self.assertEqual(nameList[index], 'secondAttr')
        index = valList.index(finalAttribute)
        self.assertEqual(nameList[index], 'finalAttribute')

    def test_insert(self):
        contents = ('<Element attr="attrValue"><Subnode>testValue</Subnode>'
                    '</Element>')
        root = Node('root')
        root2 = root.insert(contents)
        element = root2.node('Element')
        rawElement = Node.create(contents)
        self.assertEqual(element.xml(), rawElement.xml())

    def test_create(self):
        test_file_element = os.path.join(ntd.test_data_path,
                                         'some_xml_file.xml')
        fileElement = Node.create(test_file_element)
        with open(test_file_element, 'r') as myfile:
            contents = myfile.read().replace('\n', '')
        root = Node('root')
        root = root.insert(contents)
        rawElement = root.children[0]
        self.assertEqual(fileElement.xml(), rawElement.xml())

    def test_delete_attribute(self):
        tag = 'Root'
        value = '   Value   '
        anAttr = 'elValue'
        node = Node(tag, value=value, anAttr=anAttr)
        self.assertIn('anAttr', node.attributes)
        node.delAttribute('anAttr')
        self.assertNotIn('anAttr', node.attributes)

    def test_add_node(self):
        rootTag = 'Root'
        root = Node(rootTag)
        firstLevelTag = 'FirstLevel'
        firstLevel = Node(firstLevelTag)
        root += firstLevel
        self.assertIn(firstLevel, root.children)

    def test_add_nodes(self):
        rootTag = 'Root'
        root = Node(rootTag)
        firstLevelTag = 'FirstLevel'
        firstLevel = Node(firstLevelTag)
        root += firstLevel
        firstLevel2 = Node(firstLevelTag)
        root += firstLevel2
        firstLevel2ndTag = 'FirstLevel2ndTag'
        firstLevel3 = Node(firstLevel2ndTag)
        root = root + firstLevel3
        self.assertIn(firstLevel, root.children)
        self.assertIn(firstLevel2, root.children)
        self.assertIn(firstLevel3, root.children)

    def test_xml(self):
        rootTag = 'Root'
        root = Node(rootTag)
        firstLevelTag = 'FirstLevel'
        firstLevel = Node(firstLevelTag)
        root += firstLevel
        firstLevel2 = Node(firstLevelTag)
        root += firstLevel2
        firstLevel2ndTag = 'FirstLevel2ndTag'
        firstLevel3 = Node(firstLevel2ndTag)
        root += firstLevel3
        self.assertEqual(root.xml(),
                         ('<Root>\n'
                          '  <FirstLevel/>\n'
                          '  <FirstLevel/>\n'
                          '  <FirstLevel2ndTag/>\n'
                          '</Root>\n'),)

    def test_replace_node(self):
        rootTag = 'Root'
        root = Node(rootTag)
        firstLevelTag = 'FirstLevel'
        firstLevel = Node(firstLevelTag)
        root += firstLevel
        firstLevel2 = Node(firstLevelTag)
        root += firstLevel2
        firstLevel2ndTag = 'FirstLevel2ndTag'
        firstLevel3 = Node(firstLevel2ndTag)
        root.replaceNode(firstLevelTag, 1, firstLevel3)
        self.assertIn(firstLevel, root.children)
        self.assertNotIn(firstLevel2, root.children)
        self.assertIn(firstLevel3, root.children)
        self.assertEqual(len(root.children), 2)

    def test_search_node(self):
        rootTag = 'Root'
        root = Node(rootTag)
        firstLevelTag = 'FirstLevel'
        firstLevel = Node(firstLevelTag)
        root += firstLevel
        firstLevel2 = Node(firstLevelTag)
        root += firstLevel2
        firstLevel2ndTag = 'FirstLevel2ndTag'
        firstLevel3 = Node(firstLevel2ndTag)
        root += firstLevel3
        self.assertEqual(root.node(firstLevelTag,0), firstLevel)
        self.assertEqual(root.node(firstLevelTag,1), firstLevel2)

    def test_str(self):
        tag = 'Root'
        value = 'Value'
        node = Node(tag, value=value)
        self.assertEqual(str(node), '%s\n    value: [%s]' % (tag, value))

if __name__ == "__main__":
    unittest.main()
