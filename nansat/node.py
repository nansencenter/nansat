# Name:    node.py
# Purpose: Rapidly assemble XML using minimal coding.
# Authors:      Bruce Eckel, (c)2006 MindView Inc. www.MindView.net
# Contributors: Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2013
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
from __future__ import unicode_literals

import os
import re
import xml.dom.minidom as xdm


class Node(object):
    """
    Rapidly assemble XML using minimal coding.

    By Bruce Eckel, (c)2006 MindView Inc. www.MindView.net
    Permission is granted to use or modify without payment as
    long as this copyright notice is retained.

    Everything is a Node, and each Node can either have a value
    or subnodes. Subnodes can be appended to Nodes using '+=',
    and a group of Nodes can be strung together using '+'.

    Create a node containing a value by saying
    Node('tag', 'value')
    You can also give attributes to the node in the constructor:
    Node('tag', 'value', attr1 = 'attr1', attr2 = 'attr2')
    or without a value:
    Node('tag', attr1 = 'attr1', attr2 = 'attr2')

    To produce xml from a finished Node n, say n.xml() (for
    nicely formatted output) or n.rawxml().

    You can read and modify the attributes of an xml Node using
    getAttribute(), setAttribute(), or delAttribute().

    You can find the value of the first subnode with tag == 'tag'
    by saying n['tag']. If there are multiple instances of n['tag'],
    this will only find the first one, so you should use node() or
    nodeList() to narrow your search down to a Node that only has
    one instance of n['tag'] first.

    You can replace the value of the first subnode with tag == 'tag'
    by saying n['tag'] = newValue. The same issues exist as noted
    in the above paragraph.

    You can find the first node with tag == 'tag' by saying
    node('tag'). If there are multiple nodes with the same tag
    at the same level, use nodeList('tag').

    The Node class is also designed to create a kind of 'domain
    specific language' by subclassing Node to create Node types
    specific to your problem domain.

    This implementation uses xml.dom.minidom which is available
    in the standard Python 2.4 library. However, it can be
    retargeted to use other XML libraries without much effort.

    """

    def __init__(self, tag, value=None, **attributes):
        """Everything is a Node. The XML is maintained as (very efficient)
        Python objects until an XML representation is needed.

        """
        self.tag = tag.strip()
        self.attributes = attributes
        self.children = []
        self.value = value
        if self.value:
            self.value = self.value.strip()

    def getAttribute(self, name):
        """ Read XML attribute of this node. """
        return self.attributes[name]

    def setAttribute(self, name, item):
        """ Modify XML attribute of this node. """
        self.attributes[name] = item

    def delAttribute(self, name):
        """ Remove XML attribute with this name. """
        del self.attributes[name]

    def replaceAttribute(self, name, value):
        """ replace XML arrtibute of this node. """
        del self.attributes[name]
        self.attributes[name] = value

    def node(self, tag, elemNum=0):
        """ Recursively find the first subnode with this tag.

        Parameters
        ----------
        elemNum : int
            if there are several same tag, specify which element to take.

        """
        if self.tag == tag:
            return self
        ielm = 0
        for child in self.children:
            result = child.node(tag)
            if result and ielm == elemNum:
                return result
            elif result:
                ielm += 1
        return False

    def replaceNode(self, tag, elemNum=0, newNode=None):
        ''' Find the first subnode with this tag and replace with given node.

        Parameters
        ----------
        tag : str
            node tag
        elemNum : int
            number of subnode among other subnodes with similar tag

        '''
        status = False
        elemi = 0
        for i in range(len(self.children)):
            if str(self.children[i].tag) == tag:
                if elemi == elemNum:
                    self.children[i] = newNode
                    status = True
                    break
                else:
                    elemi += 1
        return status

    def delNode(self, tag, options=None):
        '''
        Recursively find nodes containing subnodes with this tag and remove
        subnodes

        Parameters
        ----------
        options : dictionary
            if there are several tags, specify a node by their attributes.

        '''
        # indeces of children to be removed
        ideleted = []

        for i, child in enumerate(self.children):
            if str(child.tag) != str(tag):
                # if child has another tag : delete children
                child.delNode(tag, options)
            elif options is None:
                # if child has this tag and no options : mark for deletion
                ideleted.append(i)
            else:
                # if child has this tag
                # and options match attributes : mark for deletion
                for key in options.keys():
                    if str(child.attributes.get(key, '')) == str(options[key]):
                        ideleted.append(i)

        # delete marked children
        for i in sorted(ideleted, reverse=True):
            self.children.pop(i)

    def nodeList(self, tag):
        '''
        Produce a list of subnodes with the same tag.

        Note
        ----
        It only makes sense to do this for the immediate
        children of a node. If you went another level down,
        the results would be ambiguous, so the user must
        choose the node to iterate over.

        '''
        return [n for n in self.children if n.tag == tag]

    def tagList(self):
        ''' Produce a list of all tags of the immediate children '''
        taglist = []
        for iTag in range(len(self.children)):
            taglist.append(str(self.children[iTag].tag))
        return taglist

    def replaceTag(self, oldTag, newTag):
        ''' Replace tag name '''
        tags = self.tagList()
        for itag, itagName in enumerate(tags):
            if itagName == oldTag:
                self.node(self.tag).children[itag].tag = newTag

    def getAttributeList(self):
        ''' get attributes and valuse from the node and return their lists '''
        nameList = []
        valList = []
        for key, val in self.attributes.items():
            nameList.append(key)
            valList.append(val)
        return nameList, valList

    def insert(self, contents):
        ''' return Node of the node with inserted <contents>'''
        dom2 = xdm.parseString(contents)
        dom1 = xdm.parseString(self.dom().toxml())
        dom1.childNodes[0].appendChild(dom1.importNode(dom2.childNodes[0],
                                                       True))
        return Node.create(dom1)

    def __getitem__(self, tag):
        '''
        Produce the value of a single subnode using operator[].
        Recursively find the first subnode with this tag.
        If you want duplicate subnodes with this tag, use
        nodeList().

        '''
        subnode = self.node(tag)
        if not subnode:
            raise KeyError
        return subnode.value

    def __setitem__(self, tag, newValue):
        '''
        Replace the value of the first subnode containing 'tag'
        with a new value, using operator[].

        '''
        assert isinstance(newValue, str), ('Value %s must be a string'
                                           % str(newValue))
        subnode = self.node(tag)
        if not subnode:
            raise KeyError
        subnode.value = newValue

    def __iadd__(self, other):
        ''' Add child nodes using operator += '''
        assert isinstance(other, Node), 'Tried to += ' + str(other)
        self.children.append(other)
        return self

    def __add__(self, other):
        ''' Allow operator + to combine children '''
        return self.__iadd__(other)

    def __str__(self):
        ''' Display this object (for debugging) '''
        result = self.tag + '\n'
        for k, v in self.attributes.items():
            result += '    attribute: %s = %s\n' % (k, v)
        if self.value:
            result += '    value: [%s]' % self.value
        return result

    # The following are the only methods that rely on the underlying
    # Implementation, and thus the only methods that need to change
    # in order to retarget to a different underlying implementation.

    # A static dom implementation object, used to create elements:
    doc = xdm.getDOMImplementation().createDocument(None, None, None)

    def dom(self):
        '''
        Lazily create a minidom from the information stored
        in this Node object.

        '''
        element = Node.doc.createElement(self.tag)
        for key, val in self.attributes.items():
            element.setAttribute(key, val)
        if self.value:
            assert not self.children, ('cannot have value and children: %s'
                                       % str(self))
            element.appendChild(Node.doc.createTextNode(self.value))
        else:
            for child in self.children:
                element.appendChild(child.dom())  # Generate children as well
        return element

    def xml(self, separator='  '):
        return self.dom().toprettyxml(separator)

    def rawxml(self):
        return str(self.dom().toxml())

    @staticmethod
    def create(dom):
        """
        Create a Node representation, given either
        a string representation of an XML doc, or a dom.

        """
        if isinstance(dom, str):
            if os.path.exists(dom):
                # parse input file
                dom = xdm.parse(dom)
            else:
                # Strip all extraneous whitespace so that
                # text input is handled consistently:
                dom = re.sub('\s+', ' ', dom)
                dom = dom.replace('> ', '>')
                dom = dom.replace(' <', '<')
                return Node.create(xdm.parseString(str(dom)))

        # To pass test for python3, decoding of bytes object is requested
        if dom.nodeType == dom.DOCUMENT_NODE:
            return Node.create(dom.childNodes[0])
        if dom.nodeName == '#text':
            return
        node = Node(dom.nodeName)
        if dom.attributes:
            for name, val in dom.attributes.items():
                node.setAttribute(name, val)
        for n in dom.childNodes:
            if n.nodeType == n.TEXT_NODE and n.wholeText.strip():
                node.value = n.wholeText
            else:
                subnode = Node.create(n)
                if subnode:
                    node += subnode
        return node
