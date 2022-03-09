#!/usr/bin/env python
import re, argparse
import itertools
import cairo 
import math
from IPython.display import display, SVG, Image
import numpy as np
import Bioinfo
import pprint
import unicodedata

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta_file", help="Input Sequence File", required=True)
	parser.add_argument("-m", "--motif_file", help="Input Motif list", required=True)
	return parser.parse_args()

args = get_args() 
m = args.motif_file
f = args.fasta_file

class Seq: 

	#Define Properties 
	def __init__ (self, m_file, seq_file):
		self.motif_file = m_file
		self.seq_file = seq_file


	def build_motif_list(self):
		''' Reads a motif file which consists one motif per line, expands each motif file into a regex char class and stores it in a list'''
		motif_list = [] #create an empty list to include all the motifs from the motif file
		#motif dict contains all the ambiguous characters as keys and the actual nucleotide replacements as values
		amb_motif_dict = {"R":["A","G"], "Y": ["C", "U", "T"] ,"S": ["C", "G"] , "W": ["A", "U"] , "K": ["G", "U"] , "M":["A", "C"] , "B": ["C", "G", "U"] , "D": ["A", "G", "U"] , "H": ["A", "C", "U"] , "V": ["A", "C", "G"] , "N": ["A", "C", "G"],"y":["c","u","t"]}
		with open (self.motif_file, "r") as m:
			char_count = 0 #maintain a count for every nucleotide
			for line in m:
				line = line.strip()
				motif_string = ""
				for char in line:
					#char_count += 1
					if char in amb_motif_dict:
						#If this is part of amb_dict, then encode the regex character class [..] here
						char_class = "["
						for element in amb_motif_dict[char]: #looping through list of values for each key in the dict
							char_class = char_class + element #Adding each element in the amb_dict to the char class
						char_class = char_class + "]"
						motif_string = motif_string + char_class
					else:
						#This is just a normal char, so append as it is..
						motif_string = motif_string + char
						#Append this motif_string and its length...this will enable us to use this string directly as a regex pattern
				char_count += 1
				motif_list.append(( '(?=' + motif_string + ')', char_count))
		print(motif_list)
				#create a new attribute calld motif_list in the object
		self.motif_list = motif_list

	
	def build_seq_list(self): 
		''' Reads each line from the sequence file and converts multi line fasta sequence to single line fasta sequence and stores the sequences in a list '''
		seq_list = Bioinfo.multi_to_single_line(f)
		self.seq_list = seq_list
		#pprint.pprint(seq_list)
		#seq_list.append(header)
		#seq_list.append(seq)

	def parse_seqfile(self):
		''' Takes the sequence list and motif list as inputs and finds motifs, replaces introns with "_" and exons with "-" in the sequence from the sequence file '''
		line_count = 0 
		header_list = []
		final_list = []
		plist= []
		seq_num = 0
		
		motif_dict = {} # key = motif, value = start and end positions of motif, 
		head = ""
		#First replace all motifs patterns in the line with a '#'
		for line in self.seq_list:
			if line.startswith(">"):
				head = line
				header_list.append(line)
				continue
			else:
				seq_num += 1
				for (motif_string, motif_identifier) in self.motif_list:
					for m in re.finditer(motif_string, line):
						plist.append((m.start(), motif_identifier, seq_num))
				#replace lower chars with a '_'
				line = re.sub(r'[a-z]', '_', line)
				#replace UPPER chars with a '-'
				line = re.sub(r'[A-Z]', '-', line)

				length = 1
				prev_char = line[0]
				start_pos = 0
				curr_char = ''
				#This will contain the final list required for drawing the figure
				
				#Iterate through the seq and mark the positions of different groups..
				for i in range(1, len(line)):
					curr_char = line[i]
					if curr_char == prev_char:
						length += 1
					else:
						final_list.append([prev_char, start_pos, length - 1])
						length = 1
						start_pos = i
						prev_char = curr_char
				final_list.append([prev_char, start_pos, length - 1])

		self.final_list = final_list
		self.plist = plist 
		import pprint
		#pprint.pprint(final_list)
		pprint.pprint(plist)
		#trial_list = [["#", 0, 50],["_", 50, 50],["-", 100 , 50],["#", 0, 30],["_", 30, 75], ["#", 105, 20], ["-", 0, 30]]
		#trial_list = [["#", 0, 30],["_", 30, 75], ["#", 105, 20]]
		#self.final_list = trial_list
		#self.plist = trial_list

    
	
class draw_fig(object): 
	''' This class is used to draw the introns, exons and motifs using the functions defined in the Seq Class. Using the object passed into this class definition, access the methods and data within my_obj (obj of Seq class) from within the draw_fig class. '''

	#Constructor
	def __init__ (self, my_obj):
		
		WIDTH, HEIGHT = 1000,1000
		self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
		self.context = cairo.Context(self.surface)

	def horizontal(self):
		'''This function is used to draw the introns'''
		
		#print("I am horizontal", self.curr_pos,self.l) 
		self.context.set_line_width(6)
		self.context.rectangle(self.curr_pos, self.y_cord + 25, self.l, 0)
		self.context.stroke()
		self.surface.write_to_png('Figure_1.png')

	def vertical(self):
		''' This function is used to draw the motifs '''
		#print("I am vertical", self.curr_pos,self.l) 
		self.context.set_line_width(2*self.l)
		self.context.stroke()
		self.context.set_source_rgb(1, 0, 0)
	
		self.context.rectangle(self.curr_pos, self.y_cord, self.l, 50)
		self.context.fill()
		self.surface.write_to_png('Figure_1.png')
		
	def box(self):
		''' This function is used to draw the exons ''' 
		#print("I am a rectangle", self.curr_pos,self.l)
		self.context.set_source_rgb(0,0,1)
		#print(self.curr_pos,self.final_pos)
		self.context.rectangle(self.curr_pos, self.y_cord, self.l, 50)
			#ctx.rectangle(x, y, width, height)
		self.context.set_line_width(6.0)
		self.context.set_source_rgb(1, 0, 0)
		self.context.stroke()

		self.surface.write_to_png('Figure_1.png')
		

	def st_line(self):
		''' This function draws motifs '''
		print("I am st_line fn:")
		print(self.pos,self.y_axis)
		self.context.rectangle(self.pos, self.y_axis, 20, 50)
		
		self.context.set_line_width(6.0)
		self.context.stroke()
		self.surface.write_to_png('Figure_1.png')

	
	

	def core_logic(self):
		''' This function is used to call the exon, intron and motif functions when required'''
		print("I am the core function")
		
		self.context.set_source_rgb(1, 0, 0)
		self.context.set_line_width(6.0)
		self.context.fill()
		curr_element = ""
		curr_pos = 0
		final_pos = 0
		y_cord = 0
		y_axis = 0
		num_count = 0


		for i in my_obj.final_list:

			curr_element = i[0]
			curr_pos = i[1]
			l = i[2]
			final_pos = l + curr_pos 
			
			self.curr_element = curr_element 
			self.curr_pos = curr_pos
			self.final_pos = final_pos 
			self.l = l
			
            #print(self.curr_element,self.)

			if self.curr_pos == 0:
				#print("here")
				y_cord = y_cord + 70
				self.context.move_to(self.curr_pos, y_cord)
				#print(self.curr_pos, y_cord)
				self.y_cord = y_cord
				self.y_axis = y_axis
			if self.curr_element == '_':
				#print("This is an intron")
				self.horizontal()
			elif self.curr_element == '-':
				#print("This is an exon")
				self.box()
			else: 
				#print ("This is a motif")
				self.vertical() 

		
		#resetting cursor to (0,0)
		my_obj2.context.move_to(0,self.y_axis)
		
		#Iterate through each element of the plist and if m_id=1, then draw a line of color x and if m_id =2, then draw a line of color y, etc..
		for index, tuple in enumerate(my_obj.plist): 
			pos =  tuple[0]
			m_id = tuple[1]
			snum = tuple[2]
			
			self.pos = pos
			self.m_id = m_id
			self.snum = snum

			#print(self.snum)	
			
			if snum == 1:
				self.y_axis = 70
				#print(self.pos,y_cord)
				if m_id == 1:
					self.context.set_source_rgb(0.5, 0, 0) 
					self.st_line()
				elif m_id == 2: 
					self.context.set_source_rgb(0, 0.5, 0) 
					self.st_line() 
				elif m_id == 3: 
					self.context.set_source_rgb(0, 0, 0.5) 
					self.st_line()
				elif m_id == 4: 
					self.context.set_source_rgb(0.5, 0.5, 0.5) 
					self.st_line()

			elif snum == 2: 
				
				self.y_axis = 140
				if m_id == 1:
					self.context.set_source_rgb(0.5, 0, 0) 
					self.st_line()
				elif m_id == 2: 
					self.context.set_source_rgb(0, 0.5, 0) 
					self.st_line() 
				elif m_id == 3: 
					self.context.set_source_rgb(0, 0, 0.5) 
					self.st_line()
				elif m_id == 4: 
					self.context.set_source_rgb(0.5, 0.5, 0.5) 
					self.st_line()

			elif snum == 3: 
				#print("I am 3rd seq")
				num_count = num_count + 1
				self.y_axis = 210
				
				if num_count == 1:
					#self.y_axis = self.y_axis + 70
					print("I am the first motif in 3rd seq ")
					if m_id == 1:
						self.context.set_source_rgb(0.5, 0, 0) 
						self.st_line()
					elif m_id == 2: 
						self.context.set_source_rgb(0, 0.5, 0) 
						self.st_line() 
					elif m_id == 3: 
						self.context.set_source_rgb(0, 0, 0.5) 
						self.st_line()
					elif m_id == 4: 
						self.context.set_source_rgb(0.5, 0.5, 0.5) 
						self.st_line()
				else: 
					#print("My count is:", num_count)
					if m_id == 1:
						self.context.set_source_rgb(0.5, 0, 0) 
						self.st_line()
					elif m_id == 2: 
						self.context.set_source_rgb(0, 0.5, 0) 
						self.st_line() 
					elif m_id == 3: 
						self.context.set_source_rgb(0, 0, 0.5) 
						self.st_line()
					elif m_id == 4: 
						self.context.set_source_rgb(0.5, 0.5, 0.5) 
						self.st_line()

			elif snum == 4: 

				self.y_axis = 280
				print("I am the first motif in 4th seq ")
				if m_id == 1:
					self.context.set_source_rgb(0.5, 0, 0) 
					self.st_line()
				elif m_id == 2: 
					self.context.set_source_rgb(0, 0.5, 0) 
					self.st_line() 
				elif m_id == 3: 
					self.context.set_source_rgb(0, 0, 0.5) 
					self.st_line()
				elif m_id == 4: 
					self.context.set_source_rgb(0.5, 0.5, 0.5) 
					self.st_line()
				
		

#Instantiating objects:
my_obj = Seq(m, f)
my_obj.build_motif_list()
my_obj.build_seq_list()
my_obj.parse_seqfile()
my_obj2 = draw_fig(my_obj)
my_obj2.core_logic()


