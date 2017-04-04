#!/usr/bin/python

# gb2tab - comprehensive GenBank format parser/extractor
#
# Copyright (C) 2004 - 2008  Rasmus Wernersson, raz@cbs.dtu.dk
#
# gb2tab is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# gb2tab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

"""
gb2tab v 1.2.1 (command line program behind the FeatureExtract webserver)

NAME
	gb2tab - extract sequence and annotation (intron/exon etc)
	         from GenBank format files.

SYNOPSIS
	gb2tab [-f 'CDS,mRNA,...'] [options...] [files...]

DESCRIPTION
	gb2tab is a tool for extracting sequence and annotation
	(such as intron / exon structure) information from GenBank
	format files.

	This tool handles overlapping genes gracefully.

	If no files are specified input is assumed to be on STDIN.
	Several GenBank files can be concatenated to STDIN.

	The extracted sequences are streamed to STDOUT with one
	entry per line in the following format (tab separated):

	name	seq	ann	com

	name:	The sequence id. See the --genename, --locustag and
		--entryname options below.

	seq:	The DNA sequence it self. UPPERCASE is used for the
		main sequence, lowercase is used for flanks (if any).

	ann:	Single letter sequence annotation. Position for position
		the annotation descripes the DNA sequence: The first
		letter in the annotation, descriped the annotation for
		the first position in the DNA sequence and so forth.

		The annotation code is defined as follows:

		FEATURE BLOCKS (AKA. "EXON BLOCKS")

		(	First position
		E	Exon
		T	tRNA exonic region
		R	rRNA / generic RNA exonic region
		P	Promotor
		X	Unknown feature type
		)	Last position

		?	Ambiguous first or last position

		[	First UTR region position
		3	3'UTR
		5	5'UTR
		]	Last UTR region position

			See also the --block-chars option, for a further
			explanation of feature blocks and exonic regions.

		INTRONS and FRAMESHIFTS

		D	First intron position (donor site)
		I	Intron position
		A	Last intron position (acceptor site)

		<	Start of frameshift
		F	Frameshift
		>	End of frameshift

		REGIONS WITHOUT FEATURES

		.	NULL annotation (no annotation).

		ONLY IN FLANKING REGIONS:

		+	Other feature defined on the SAME STRAND
			as the current entry.
		-	Other feature defined on the OPPOSITE STRAND
			relative to the current entry.
		#	Multiple or overlapping features.

		A..Z:	Feature on the SAME STRAND as the current entry.
		a..z:	Feature on the OPPOSITE STRAND as the current entry.

			See the -e option for a description of which features
			are annotated in the flanking regions.

			The options --flank_ann_full (default) and
			--flank_ann_presence determine if full annotation
			(+upper/lower case) or annotation of presence/absence
			(+/- and #) is used.


	com:	Comments (free text). All text, extra information etc
		defined in the GenBank files are concatenated into a single
		comment.

		The following extra information is added by this program:

		*) GenBank accession ID.
		*) Source (organism)
		*) Feature type (e.g. "CDS" or "rRNA")
		*) Strand ("+" or "-").
		*) Spliced DNA sequence. Simply the DNA sequence defined
		   by the JOIN statement.
		   This is provied for two reasons. 1) To overcome negative
		   frameshifts. 2) As an easy way of extracting the sequence
		   of the spliced producted. See also the --splic_always and
		   --flank_splic options below.
		*) Spliced DNA annotation.

OPTIONS
	The following options are available.

	-f X, --feature_type=X
		Define which feature type(s) to extract.

		Default is 'CDS' which is the most general way
		to annotate protein coding genes.

		Multiple features can be selected by specifying a comma
		separated list - for example "CDS,rRNA,tRNA".

		Special keywords:

		ALL: 	Using the keyword "ALL", will extend the list to all
			feature types listed in each GenBank file.
			Please notice: This can occationally lead to problems
			in files that use multiple feature types to cover the
			same actual feature (e.g uses both "gene" and "CDS").

		MOST:	Covers the following feature types:

			CDS,3'UTR,5'UTR,
			promoter,-35_signal,-10_signal,RBS,
			rRNA,tRNA,snoRNA,scRNA,misc_RNA,
			misc_feature

		The keyword can be also be included in the user specified list.

		For example "MOST,novel_feature" will construct a list containing
		the list mention above + the new feature type "novel_feature".

	-e X, --flank_features=X
		Define which features to annotate in flanking regions.

		The scheme for specifying features is the same as in the
		-f option (see above).

		The default value is "MOST".

		If no flanking regions are requested (see options -b and -a
		below) this option is ignored.

	-i, --intergenic
		Extract intergenic regions. When this options is used all
		regions in between the features defined with the -f options
		in extracted rahter than the features themselves.

		Please notice that features specified using the -e options
		may be present in the intergenic regions.

		Intergenic regions will always be extracted from the "+" strand.

	-s, --splice
		For intron containing sequences output the spliced version as
		the main result (normally this information goes into the
		comments). If this options is used the full length product will
		be added to the comments instead.

		Using this option will force the inclusion of flanks (if any)
		in the spliced product. See also option --flank_splic.

	-x, --spliced_only
		Only output intron containing sequences. Can the used in
		combination with the -s option.

	-b X, --flank_before=X
		Extract X basepairs upstream of each sequence.

	-a X, --flank_after=X
		Extract X basepairs downstream of each sequence.

	-h, --help
		Print this help page and exit.

	-n, --dry-run
		Run through all extraction steps but do not output any
		data. Useful for debugging bad GenBank files in combination
		with the verbose options.

	-v, --verbose
		Output messages about progess, details about the GenBank
		file etc. to STDERR. Useful for finding errors.

	-q, --quiet
		Suppress all warnings, error messages and verbose info.
		The exit value will still be non-zero if an error is
		encountered.

	--flank_ann_presence
		Annotate presence/absence and relative strandness of
		features in the flanking regions.

		Features - of any kind - are annotated with "+" if they are
		on the SAME STRAND as the extratced feature, and "-" if they
		are on the OPPOSITE STRAND. "#" marks regions covered by
		multiple features.

		This option is very useful for use with OligoWiz-2.0
		(www.cbs.dtu.dk/services/OligoWiz2).

	--flank_ann_full
		Default: Include full-featured annotation in the flanking regions.

		Features on the SAME STRAND as the extracted is uppercase -
		features on the OPPOSITE STRAND is lowercase.

		In case of regions covered by multiple features, the
		feature defined FIRST by the -e option has preference.

	--flank_splic
		Also include flanking regions in the spliced product.
		Default is to ignore flanks.

	--splic_always
		Include spliced producted for ALL entries.
		Default is to only print spliced product information for
		intron/frameshift containing entries.

	--frameshift=X
		"Introns" shorter than X bp (default 15bp) are considered
		frameshifts. This includes negative frameshifts.

	--block-chars=XYZ|"Feat1=XYZ,Feat2=ZYX,..."
		Specify which characters to use for annotation of the
		extracted feature types. For spliced feature (e.g CDS)
		each exonic block is annotated using the specified characters.

		Three characters must be supplied (for each feature type):
		First position, internal positions, last position.

		For example the string "(E)" will cause a 10bp feature block
		(e.i a CDS exon block) to be annotated like this: (EEEEEEEE)

		Introns are filled in as DII..IIA

		By default the program determine the annotation chars to be
		based on the type of feature being extracted:

		(E)	CDS, mRNA
		(T)	tRNA
		(R)	rRNA, snoRNA, snRNA, misc_RNA, scRNA
		(P)	promotor
		[5]	5'UTR
		[3]	3'UTR

		(X)	Everything else.

		This table can be expanded (and overwritten) by supplying a
		list of relations between feature type ans block chars.
		For example:

		--block-chars="mRNA=[M],gene=<G>,repeat=QQQ"

	--genename
		Try to extract the gene name from the /gene="xxxx"
		tag (this is usually the classical gene name, e.g. HTA1)
		If this is not possible fall back to 1) locustag
		or 2) entryname (see below).

	--locustag
		Try to extract the locus tag (usually the systematic
		gene name) from the /locus_tag="xxxx" tag. Fall back
		to using the entryname if not possible (see below).

		This is the default behavior.

	--entryname
		Use the main GenBank entry name (the "LOCUS" name) as
		the base of the sequence names.

KNOWN ISSUES
	This program DOES NOT support entries which spans multiple
	GenBank files. It is very unlikely this will ever be supported.
	(Please notice that the webserver version supports expanding
	reference GenBank entries to the listed subentries automatically).

REFERENCE
	Rasmus Wernersson, 2005.
	"FeatureExtract - extraction of sequence annotation made easy".
	Nucleic Acids Research, 2005, Vol. 33, Web Server issue W567-W569

WEB
	http://www.cbs.dtu.dk/services/FeatureExtract

	The webpage contains detailed instructions and examples.
	The most recent version of this program is downloadable
	from this web address.

AUTHOR
	Rasmus Wernersson, raz@cbs.dtu.dk

	Oct-Dec 2004
	Jan-Mar 2005
	Aug     2005
	Sep     2008 - bugfix + better IUPAC support

"""


import sys,string,re

# The Exception is used for communicating the end of the input stream has be reached
class ExEndOfFile(Exception):
	pass

# The class is basically just a slightly intelligent record (it knows how to be
# sorted) that stores the processed information for each feature block.
class Rec:
	def __init__(self):
		self.first      = -1
		self.last       = -1
		self.comp       = False
		self.strand     = "+"
		self.strSeq     = ""
		self.strAnn     = ""
		self.strSpliced = ""
		self.strSpcAnn  = ""
		self.featType   = ""
		self.comment    = ""
	def __cmp__(self,other):
		return cmp(self.first,other.first)
	def __str__(self):
		return "first: %i, last: %i" % (self.first,self.last)

# The MOST keyword table from the man page can be pasted in here:
KeyWord_MOST = """
			CDS,3'UTR,5'UTR,
			promoter,-35_signal,-10_signal,RBS,
			rRNA,tRNA,snoRNA,scRNA,misc_RNA,
			misc_feature
""".replace("\n","").replace("\t","").strip()

# The exon-char table from the man page can be inserted here
exon_chr_table = """
		(E)	CDS, mRNA
		(T)	tRNA
		(R)	rRNA, snoRNA, snRNA, misc_RNA, scRNA
		(P)	promoter
		[5]	5'UTR
		[3]	3'UTR

"""

d_echr = {"":"(X)"}
for line in exon_chr_table.split("\n"):
	line = line.strip()
	if not line: continue
	tokens = line.split()
	echr = tokens[0]
	for key in tokens[1:]:
		key = key.strip(",")
		d_echr[key] = echr

# Reverse compliment DNA.
#
# IUPAC alphabet:
# --------------
# Code    Meaning             Compliment
# A	  A                   T
# C       C                   G
# G       G                   C
# T       T                   A
#
# R       A, G (puRine)       T, C (Y)
# Y       C, T (pYrimidine)   G, A (R)
# M       A, C (aMino)        T, G (K)
# K       G, T (Keto)         C, A (M)
# S       C, G (Strong)       G, C (S)
# W       A, T (Weak)         T, A (W)
#
# B       C, G, T (not A)     G, C, T (V)
# D       A, G, T (not C)     T, C, A (H)
# H       A, C, T (not G)     T, G, A (D)
# V       A, C, G (not T)     T, G, C (B)
#
# N       A, C, G, T (aNy)    T, G, C, A (N)
alphaPlus  = "ACGTRYMKSWBDHVN"
alphaMinus = "TGCAYRKMSWVHDBN"

def rcDna(seq):
	lseq = list(seq)
	lseq.reverse()
	rev = "".join(lseq)
	rcomp = rev.translate(string.maketrans( alphaPlus+alphaPlus.lower() , alphaMinus+alphaMinus.lower() ))
	return rcomp

# Reverse compliment the annotation string. Notice how block start and end markers
# like "(" and ")" are swapped around.
def rcAnn(ann):
	lann = list(ann)
	lann.reverse()
	rev = "".join(lann)
	rcomp = rev.translate(string.maketrans("{}[53](EE)BMS<F>DIA+-#","}{]35[)EE(SMB>F<AID-+#") )
#	rcomp = rev.translate(string.maketrans("[53](EC)BMS<F>DIA+-#","]35[)CE(SMB>F<AID-+#") )
	return rcomp

################
# PARSING STEP 1: Run through the lines in the actual GenBank format file/stream, and build
#                 feature block records for later use.
#
#                 Notice that this function works on an input stream (not a file) and the parsing
#                 has therefore been built to be single-pass.
################
def build_FeatureList(in_stream,opt):
	l_seq = []

	d_head = {}
	d_feat = {}

	cur_feature_type = ""
	l = []

	inHeader   = True
	inFeatures = False
	inSeq      = False

	locusLine  = ""

	while True:
		line = in_stream.readline()

		# Check for premature termination of the file
		if not line:
			in_stream.close()

			# Ternimation is premature, if we are actually in the process of parsing...
			if locusLine:
				sys.stderr.write("WARNING: Unexpected EndOfFile for GenBank entry: %s" % locusLine)

			# ... Otherwise is a perfecly normal situation.
			raise ExEndOfFile()

		if line.startswith("//"):
			# End of this GenBank entry. The file/stream may contains multiple entries.
			break

		elif inHeader:
			if line.startswith("FEATURES"):
				if opt.debug: sys.stderr.write("Entering FEATURES section\n")
				inHeader   = False
				inFeatures = True
			else:
				if line.startswith("LOCUS"):
					locusLine = line

				if not line.startswith(" "):
					name = line[:12].strip()
					val  = line[12:-1]
					d_head[name] = val
		elif inSeq:
			tokens = line.split()
			l_seq.extend(tokens[1:])

		else:
			if not inFeatures:
				# Wait for start of the ORIGIN (sequence containing) block
				if line.startswith("ORIGIN"):
					inSeq = True
				continue

			# Test if we have left the FEATURE block
			if not line.startswith(" "):
				# Directly into the ORIGIN (sequence) block?
				if line.startswith("ORIGIN"):
					inSeq = True
				# .. Or perphaps just into an unknown type of data block
				inFeatures = False
				continue

			name = line[:21].strip()
			val  = line[21:-1]

			if name == "":
				l.append(val)

			else:
				# We have reached a new feature block
				if not d_feat.has_key(name):
					d_feat[name] = []

				# The trick is here that we put the new (empty)
				# list object directly into the feature table.
				# This way lines added to "l" later will also
				# be found throught the reference to "l" in the
				# feature table.
				cur_feature_type = name
				l = []
				d_feat[cur_feature_type].append(l)
				l.append(val)

	# Parsing of this entry is complete.
	seq = "".join(l_seq)
	return (d_feat, d_head, seq)

################
# PARSING STEP 2: Run through the featurelist generated (and later filtered) from the parsing of the
#                 GenBank file. In this step the annotations strings used outputted later is constructed
#                 and, most importantly, the annotation used in flanking regions are added to the
#                 whole sequence annotation strings.
#
################
def parse_FeatureList(featureList, opt, seq, annPlus, annMinus, annCombined,curFeat,locusName):
	# Determine exon chars for this feature type
	if opt.blockChars:
		# Custom defined - overrule the look up table
		exChr = opt.blockChars
	else:
		# Default - base exon chars on the feature type
		if d_echr.has_key(curFeat):
			exChr = d_echr[curFeat]
		else:
			exChr = d_echr[""]

	ex1, ex2, ex3 = exChr[0], exChr[1], exChr[2]

	# Run through the feature list
	lstExtracted = []
	for ent in featureList:
		# Pos line is the very first lines, until a new tag is found
		pos_line = ""
		for i in range(0,len(ent)):
			line = ent[i].strip()
			if line.startswith("/"): break
			pos_line += line

		# Check for unsupported stuff in the JOIN statement
		pos_err = ""
		if ":" in pos_line:
			pos_err = "Features spanning multiple GenBank entries are not supported."
#		elif (">" in pos_line or "<" in pos_line):
#			pos_err = "Ambiguous positions - Cannot extract this feature."

		if pos_err:
			s = "SKIPPING feature (GB entry: %s) - reason: %s\nThe JOIN statement is: %s\n"
			s = s % (locusName, pos_err, pos_line)

			if not opt.quiet:
				sys.stderr.write(s)

			continue

		# Sum up everything after the position as a comment
		comment = "".join(ent[i:])

		comp = 0

		if pos_line.startswith("complement("):
			comp = 1
			pos_line = pos_line[len("complement("):-1]

		if pos_line.startswith("join("):
			pos_line = pos_line[len("join("):-1]

		pos_list = pos_line.split(",")

		if comp:
			annWork = annMinus
		else:
			annWork = annPlus

		errMsg    = ""
		l_spliced = []
		l_spc_ann = []
		pre_pos = -1
		first   = -1
		last    = -1
		for pos in pos_list:
			tokens = pos.split("..")
			try:
				ambStart = ambEnd = False

				strPos = tokens[0]
				if "<" in strPos:
					ambStart = True
					strPos = strPos.replace("<","")

				start = int(strPos) - 1 # Adjust to zero-based

				# Handle a rare situation where the first exon is only 1bp long
				# (happens in the Yeast genome - chromosome 9).
				if len(tokens) == 1:
					end = start
				else:
					strPos = tokens[1]
					if ">" in strPos:
						ambEnd = True
						strPos = strPos.replace(">","")

					end   = int(strPos) - 1 # Adjust to zero-based

				#Record first and last pos
				if (first == -1): first = start
				last = end
			except:
				# Handle a special case: a gene with partial DNA inversion (entry ECP15BG)
				if "complement" in pos:
					errMsg = "Features with partial DNA inversion not supported"
				else:
					errMsg = "Error parsing JOIN stament."
				break

			# Handle a special case: circular genes (See gb entry AACOP - first CDS)
			if start < first:
				errMsg = "Features with circular sequence segments not supported."
				break

			# Fill in intron / frameshift
			if pre_pos <> -1:
				if (start - pre_pos) < opt.frameShiftCutOff:
					chrB = "<"
					chrF = "F"
					chrE = ">"
				else:
					chrB = "D"
					chrF = "I"
					chrE = "A"

				for i in range(pre_pos+1,start):
					annWork[i] = chrF
				annWork[pre_pos+1] = chrB
				annWork[start - 1] = chrE

			# Record position of last nucleotide in the exon
			pre_pos = end

			# Fill in exon
			for i in range(start,end+1):
				annWork[i] = ex2

			# Mark exon start and end
			if ambStart: annWork[start] = "?"
			else:        annWork[start] = ex1

			if ambEnd:   annWork[end]   = "?"
			else:        annWork[end]   = ex3

			# Add exon sequence to spliced product
			l_spliced.append(seq[start:end+1])
			l_spc_ann.extend(annWork[start:end+1])

		# Skip bad sequenes
		if errMsg:
			if not opt.quiet:
				s = "SKIPPING feature (GB entry: %s) - reason: %s \nThe JOIN statement is: %s\n"
				sys.stderr.write(s % (locusName,errMsg,pos_line))
			continue

		# post-process spliced data
		strSpliced = "".join(l_spliced)
		strSpcAnn  = "".join(l_spc_ann)
		if comp:
			strSpliced = rcDna(strSpliced)
			strSpcAnn  = rcAnn(strSpcAnn)

		# Post process seq / ann data
		strAnn = "".join(annWork[first:last+1])
		strSeq = (seq[first:last+1]).upper()
		if comp:
			strCompAnn = strAnn     #May be needed for flank annotation
			strSeq = rcDna(strSeq)
			strAnn = rcAnn(strAnn)

		#Put aside data for later use
		r = Rec()

		r.first      = first
		r.last       = last
		r.comp       = comp

		r.strSeq     = strSeq
		r.strAnn     = strAnn

		r.strSpliced = strSpliced
		r.strSpcAnn  = strSpcAnn
		r.featType   = curFeat

		r.comment    = comment

		if comp:
			r.strand = "-"
		else:
			r.strand = "+"

		lstExtracted.append(r)

		if opt.flankFullAnn:
			if comp: fillAnn = strCompAnn.lower()
			else:    fillAnn = strAnn

			j = 0
			for i in range(first,last+1):
				annCombined[i] = fillAnn[j]
				j += 1
		else:
			if comp: chrF = "-"
			else:    chrF = "+"
			for i in range(first,last+1):
				c = annCombined[i]
				if c <> ".":
					annCombined[i] = "#"
				else:
					annCombined[i] = chrF

	return lstExtracted

#####################
# PARSING CONTROLLER: This function controls the actual parsing of the GenBank data (see
#                     function "build_FeatureList"), the generation of annotation (see
#                     function "parse_FeatureList") and takes care of the output of the data.
#
#####################
def parse_GBstream(in_stream,out_stream,opt):
	# Parse the actual GenBank file
	try:
		d_feat, d_head, seq = build_FeatureList(in_stream,opt)
	except ExEndOfFile:
		# End of file/stream reached? This is NOT an actual error.
		return
	try:
		locusName = d_head["LOCUS"].split()[0]
	except:
		locusName = "<UNKNOWN LOCUS>"

	try:
		accession = d_head["ACCESSION"]
	except:
		accession = "<UNKNOWN ACCESSION ID>"

	try:
		organism  = d_head["SOURCE"]
	except:
		organism  = "<UNKOWN ORGANISM>"

	# Verbose info: Breakdown of features
	if opt.verbose:
		s = "Breakdown of features in GenBank entry %s\n" % locusName
		for key in d_feat.keys():
			l = d_feat[key]
			s += "%i:\t%s\n" % (len(l),key)
		sys.stderr.write(s)

	# Test if this is a GenBank entry without sequence information
	if not seq:
		if not opt.quiet:
			sys.stderr.write("SKIPPING entry '%s'. No sequence data ('ORIGIN' block missing).\n" % locusName)
		return

	# Prepare analysis of the parsed feature entries
	annPlus     = ["."] * len(seq)
	annMinus    = ["."] * len(seq)
	annCombined = ["."] * len(seq)

	lstExtracted = []

	# Build a list of wanted feature for 1) extraction and 2) annotation in flanks
	featExt = opt.featureType.split(",")
	if "ALL" in featExt: featExt = d_feat.keys()

	featAnn = opt.featureAnn.split(",")
	if "ALL" in featAnn: featAnn = d_feat.keys()

	if opt.debug:
		s = "featExt: %s\nfeatAnn: %s\n" % (",".join(featExt),",".join(featAnn))
		sys.stderr.write(s)

	noFlanks = (opt.flankBefore == 0) and (opt.flankAfter == 0)

	# Loop through the wanted features
	for curFeat in featExt:
		if d_feat.has_key(curFeat):
			featLst = d_feat[curFeat]
		else:
			featLst = []

		lstFeatExt = parse_FeatureList(featLst,opt,seq,annPlus,annMinus,annCombined,curFeat,locusName)
		lstExtracted.extend(lstFeatExt)

	# Loop through annotate only features - but only if flanks are extracted or intergenic regions are wanted
	if (not noFlanks) or (opt.interGenic):
		for curFeat in featAnn:
			if d_feat.has_key(curFeat):
				featLst = d_feat[curFeat]
			else:
				featLst = []

			# Avoid redundant work
			if not curFeat in featExt:
				parse_FeatureList(featLst,opt,seq,annPlus,annMinus,annCombined,curFeat,locusName)

	# Test if any features where extracted at all
	if not lstExtracted:
		if not opt.quiet:
			sys.stderr.write("SKIPPING. No features of type '%s' in GenBank entry %s\n" % (opt.featureType,locusName))
		return

	# Process the list of annotated features
	numExtracted = 0
	strAnnComb = "".join(annCombined)
	seqRc = rcDna(seq)

	# Sort relative to position
	lstExtracted.sort()

	# Do we want the intergenic regions?
	if opt.interGenic:
		# Notice: The list of primary features are sorted by first positon now
		lstInterGenic = []
		prvR = None
		max_last = -1
		for r in lstExtracted:

#			print max_last,r.first, r.last, r.comment

			if prvR:
				if r.first > max_last:
					newR = Rec()
#					newR.first = prvR.last +1 # DON'T: This will give a problem with overlapping / internal genes ...
					newR.first = max_last +1  # ... DO this instead.
					newR.last = r.first -1
					newR.strSeq  = seq[newR.first:newR.last+1].upper()
					newR.strAnn  = strAnnComb[newR.first:newR.last+1]
					newR.featType = "Intergenic"
					#print newR
					lstInterGenic.append(newR)
			prvR = r
			max_last = max(max_last,prvR.last)

		lstExtracted += lstInterGenic
		lstExtracted.sort()

	# Run through the features to output
	for r in lstExtracted:
		flankBefore = opt.flankBefore
		flankAfter  = opt.flankAfter
		if r.comp:
			flankBefore, flankAfter = flankAfter, flankBefore

		startPos = max(0,r.first-flankBefore)

		# Determine sequence name according to the options and the information
		# which actually available.
		seqName = ""
		if opt.nameLocustag or opt.nameGene:
			gene     = ""
			locustag = ""

			posSt = r.comment.find('/gene="')
			if posSt > -1:
				posSt += len('/gene="')
				gene = r.comment[posSt:r.comment.find('"',posSt)]

			posSt = r.comment.find('/locus_tag="')
			if posSt > -1:
				posSt += len('/locus_tag="')
				locustag = r.comment[posSt:r.comment.find('"',posSt)]

			if opt.nameGene and gene:
				seqName = gene
			else:
				seqName = locustag

			if seqName and opt.flankBefore <> 0:
				posRel = 0 - opt.flankBefore
				seqName = "%s_%i" % (seqName, posRel)

		# Fall back to a sequence name consisting of the LOCUS name + position.
		if not seqName:
			seqName = "%s_%i" % (locusName, startPos + 1) # Re-convert to 1-based

		# Handle flanking regions
		seq5prime = seq[startPos:r.first]
		seq3prime = seq[r.last+1:r.last+flankAfter+1]

		ann5prime = strAnnComb[startPos:r.first]
		ann3prime = strAnnComb[r.last+1:r.last+flankAfter+1]

		if r.comp:
			seq5prime, seq3prime = rcDna(seq3prime), rcDna(seq5prime)
			ann5prime, ann3prime = rcAnn(ann3prime), rcAnn(ann5prime)

		# Main sequence and annotaion
		totalSeq = seq5prime + r.strSeq + seq3prime
		totalAnn = ann5prime + r.strAnn + ann3prime

		# Always include info about strand + GenBank entry ID
		comment = r.comment
		comment += ' /GenBank_acc="%s";' % (accession)
		comment += ' /Source="%s";' % organism
		comment += ' /feature_type="%s";' % (r.featType)
		comment += ' /strand="%s";' % (r.strand)

		# Spliced product + annotation
		noIntrons = (r.strSpcAnn == r.strAnn)
		if noIntrons and (not opt.splicAlways):
			pass
		else:
			# Handle flanking annotation
			if opt.flanksInSplic or opt.doSplice:
				r.strSpliced = seq5prime + r.strSpliced.upper() + seq3prime
				r.strSpcAnn  = ann5prime + r.strSpcAnn  + ann3prime

			# Should the spliced or the full product be in the comments?
			if opt.doSplice:
				comment += ' /full_product="%s"; /full_annotation="%s";' % (totalSeq,totalAnn)
			else:
				comment += ' /spliced_product="%s"; /spliced_annotation="%s";' % (r.strSpliced,r.strSpcAnn)

		# Output entire entry line
		if not opt.dryRun:

			if noIntrons and opt.splicedOnly:
				continue

			if opt.doSplice:
				#print "\t".join([seqName,r.strSpliced,r.strSpcAnn,comment])
				out_stream.write("\t".join([seqName,r.strSpliced,r.strSpcAnn,comment]))
			else:
				#print "\t".join([seqName,totalSeq,totalAnn,comment])
				out_stream.write("\t".join([seqName,r.strSpliced,r.strSpcAnn,comment]))

		numExtracted += 1

	if opt.verbose:
		s = "# sequences extratced: % i [Feature type: %s]\n" % (numExtracted, opt.featureType)
		sys.stderr.write(s)


###########################
# READ COMMAND LINE OPTIONS
###########################
def getOptions():
	# Quick hack to overrule the -h and --help feature
	# build into the optpase module
	if "-h" in sys.argv or "--help" in sys.argv:
		print __doc__
		sys.exit(0)

	from optparse import OptionParser
	parser = OptionParser()

	# Feature types to extract and annotate
	parser.add_option("-f","--feature_type",   type="string",  dest="featureType", default="CDS")
	parser.add_option("-e","--flank_features", type="string",  dest="featureAnn",  default="MOST")
	parser.add_option("-i","--intergenic", action="store_true", dest="interGenic", default=False)

	# Flanks
	parser.add_option("-b","--flank_before",type="int", dest="flankBefore", default=0)
	parser.add_option("-a","--flank_after", type="int", dest="flankAfter",  default=0)

	parser.add_option("--flank_ann_full",     action="store_true",  dest="flankFullAnn",  default=True)
	parser.add_option("--flank_ann_presence", action="store_false", dest="flankFullAnn",  default=True)
	parser.add_option("--flank_splic",        action="store_true",  dest="flanksInSplic", default=False)

	# Output options
	parser.add_option("-n","--dry-run", action="store_true", dest="dryRun",  default=False)
	parser.add_option("-v","--verbose", action="store_true", dest="verbose", default=False)
	parser.add_option("-q","--quiet",   action="store_true", dest="quiet",   default=False)

	# Exon chars
	parser.add_option("--block-chars", type="string", dest="blockChars", default="")

	# Frameshift
	parser.add_option("--frameshift", type="int", dest="frameShiftCutOff",  default=15)

	# Splicing
	parser.add_option("-s","--splice",       action="store_true", dest="doSplice",    default=False)
	parser.add_option("-x","--spliced_only", action="store_true", dest="splicedOnly", default=False)
	parser.add_option("--splic_always",      action="store_true", dest="splicAlways", default=False)

	# Naming
	parser.add_option("--locustag",  action="store_true", dest="nameLocustag", default=False)
	parser.add_option("--genename",  action="store_true", dest="nameGene",     default=False)
	parser.add_option("--entryname", action="store_true", dest="nameEntry",    default=False)

	# Debug
	parser.add_option("-d","--debug", action="store_true", dest="debug", default=False)

	(options, args) = parser.parse_args()

	# Naming - there can be only one.... option
	if options.nameEntry:
		options.nameGene     = False
		options.nameLocustag = False
	elif options.nameGene:
		options.nameLocustag = False
		options.nameEntry    = False
	else:
		#Default
		options.nameGene     = False
		options.nameEntry    = False
		options.nameLocustag = True

	# Output options
	if options.quiet:
		options.verbose = False

	# Feature block chars
	# Examples og valid input
	# 1) "(E)"
	# 2) "CDS=(E),RBS=<B>,repeat=RRR"
	blkErr = False
	if options.blockChars:
		if   len(options.blockChars) < 3:
			blkErr = True

		elif len(options.blockChars) > 3:  # Indicates multi feature string
			for ftSet in options.blockChars.split(","):
				tokens = ftSet.split("=")
				try:
					featType, blockChrs = tokens
					if len(blockChrs) != 3:
						blkErr = True
						break
					else:
						d_echr[featType] = blockChrs
				except:
					blkErr = True
					break
			# Now all informations is transfered to the look up table, and
			# the blockChars variable must be cleared to indicate this.
			if not blkErr: options.blockChars = ""
	if blkErr:
		sys.stderr.write("ERROR: Malformed --block-chars statement: '%s' \n" % options.blockChars)
		sys.exit(1)

	# Feature extraction - expand keywords ("ALL" must be handled later")
	options.featureType = options.featureType.replace("MOST",KeyWord_MOST)
	options.featureAnn  = options.featureAnn.replace("MOST",KeyWord_MOST)

	return options, args

if __name__ == "__main__":
	# Echo a message to STDERR if invoked with no parameters
	if len(sys.argv) == 1:
		sys.stderr.write("Reading from STDIN. Use -h for help.\n")

	# Parse options
	opt, filelist = getOptions()

	# Run through a list of files - if any
	for fn in filelist:
		try:
			in_stream = file(fn,"rU")
			if opt.verbose:
				sys.stderr.write('Reading file "%s"\n' % fn)

			# May contain multiple entries - keep the file open until all is read
			while not in_stream.closed:
				parse_GBstream(in_stream,sys.stdout,opt)

		except Exception, e:
			s = "Error parsing file: %s\n" % fn
			s += "Error message: \n%s\n" % e

			if not opt.quiet:
				sys.stderr.write(s)

			sys.exit(1)

	if len(filelist) > 0: sys.exit(0) # All done

	# No files - Read from STDIN
	while not sys.stdin.closed:
		try:
			parse_GBstream(sys.stdin,sys.stdout,opt)
		except Exception, e:
			s ="Error parsing steam from STDIN.\n"
			s += "The problem appears to be located before the following lines:\n"

			if not sys.stdin.closed:
				s += "-----> start of citation\n"
				for i in range(1,10):
					s += sys.stdin.readline()
				s += "<----- end of citation\n"

			s += "Error message: \n%s\n" % e

			if not opt.quiet:
				sys.stderr.write(s)

			sys.exit(1)
