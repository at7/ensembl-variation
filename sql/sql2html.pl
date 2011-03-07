#!/bin/perl
# 1st Feb 2011
# Generate an HTML documentation page from an SQL file.
#
# It needs to have a "javascript like" documentation above each table. e.g.:
####################################################################################
#/**
#@table variation

#@desc This is the schema's generic representation of a variation.

#@column variation_id				Primary key, internal identifier.
#@column source_id					Foreign key references to the @link source table.
#@column name								Name of the variation. e.g. "rs1333049".
#@column validation_status	Variant discovery method and validation from dbSNP.
#@column ancestral_allele		Taken from dbSNP to show ancestral allele for the variation.
#@column flipped						This is set to 1 if the variant is flipped.
#@column class_so_id				Class of the variation, based on the Sequence Ontology.

#@see variation_synonym
#@see flanking_sequence
#@see failed_variation
#@see variation_feature
#@see variation_group_variation
#@see allele
#@see allele_group_allele
#@see individual_genotype_multiple_bp
#*/
#
#
#create table variation (
#		variation_id int(10) unsigned not null auto_increment, # PK
#		source_id int(10) unsigned not null, 
#		name varchar(255),
#		validation_status SET('cluster','freq','submitter','doublehit','hapmap','1000Genome','failed','precious'),
#		ancestral_allele text,
#		flipped tinyint(1) unsigned NULL DEFAULT NULL,
#		class_so_id ENUM('SO:0001483','SO:1000002','SO:0000667') DEFAULT 'SO:0001059', # default to sequence_alteration
#
#		primary key( variation_id ),
#		unique ( name ),
#		key source_idx (source_id)
#);
########################################################################################
#
# /** and */ : begin and end of the document block
# @set       : tag to create a set of tables
# @table     : name of the sql table
# @desc      : description of the role/content of the table, set or info tags
# @column    : column_name [tab(s)] Column description. Note: 1 ligne = 1 column
# @see       : tables names linked to the described table
# @link      : Internal link to an other table description. The format is ... @link table_name ...
# @info      : tag to describe additional information about a table or a set of tables


use strict;
use POSIX;
use Getopt::Long;

###############
### Options ###
###############
my ($sql_file,$html_file,$set_flag,$sort_tables,$help);

usage() if (!scalar(@ARGV));
 
GetOptions(
    'i=s' => \$sql_file,
    'o=s' => \$html_file,
    'show_set=i' => \$set_flag,
    'sort=i' => \$sort_tables,
		'help!' => \$help
);

usage() if ($help);

if (!$sql_file) {
	print "Error! Please give a sql file using the option '-i' \n";
	exit;
}
if (!$html_file) {
	print "Error! Please give an output file using the option '-o'\n";
	exit;
}

#$set_flag ||= 1;
$set_flag = 1 if (!defined($set_flag));
$sort_tables = 1 if (!defined($sort_tables));


##############
### Header ###
##############
my $header = qq{
<html>
<head>
<meta http-equiv="CONTENT-TYPE" content="text/html; charset=utf-8" />
<title>e! Variation schema </title>

<script language="Javascript" type="text/javascript">
	// Function to show/hide the columns table
	function show_hide (param) {
		div   = document.getElementById('div_'+param);
		alink = document.getElementById('a_'+param);
		if (div.style.display=='inline') {
			div.style.display='none';
			alink.innerHTML='Show';
		}
		else {
			if (div.style.display=='none') {
				div.style.display='inline';
				alink.innerHTML='Hide';
			}
		}
	}
</script>

</head>

<body>
<h1>Ensembl Variation Schema Documentation</h1>

<h2>Introduction</h2>

<p>
This document gives a high-level description of the tables that
make up the Ensembl variation schema. Tables are grouped into logical
groups, and the purpose of each table is explained. It is intended to
allow people to familiarise themselves with the schema when
encountering it for the first time, or when they need to use some
tables that they have not used before. Note that this document
makes no attempt to enumerate all of the names, types and contents of
every single table.
</p>

<p>
This document refers to version <strong>61</strong> of the Ensembl
variation schema. 
</p>

<p>
A PDF document of the schema is available <a href="variation-database-schema.pdf">here</a>.
</p>
};


##############
### Footer  ##
##############
my $footer = qq{
</body>
</html>};



################
### Settings  ##
################
my %display_col = ('Show' => 'none', 'Hide' => 'inline');
my $documentation = {};
my $tables_names = {'default' => []};

my $in_doc = 0;
my $in_table = 0;
my $set = 'default';
my $table = '';
my $info = '';
my $nb_by_col = 15;
my $count_sql_col = 0;
my $tag_content = '';
my $tag = '';
my $display = 'Show';


#############
## Parser  ##
#############

open SQLFILE, "< $sql_file" or die "Can't open $sql_file : $!";
while (<SQLFILE>) {
	chomp $_;
	next if ($_ eq '');
	
	# Verifications
	if ($_ =~ /^\/\*\*/)  { $in_doc=1; next; }  # start of a table documentation
	if ($_ =~ /^\s*create\stable\s(if\snot\sexists\s)?(\w+)/i) { # start to parse the content of the table
		if ($2 eq $table) { 
			$in_table=1; 
		}
		else { 
			print STDERR "The documentation of the table $2 has not be found!\n";
		}
		next;
	}	
	next if ($in_doc==0 and $in_table==0);
	
	my $doc = $_;
	
	## Parsing of the documentation ##
	if ($in_doc==1) {
		# Set of tables name
		if ($doc =~ /^\@set\s*(\w+)/i and $set_flag == 1) {
			$set = $1;
			$tables_names->{$set} = [];
			next;
		}		
		# Table name
		elsif ($doc =~ /^\@table\s*(\w+)/i) {
			$table = $1;
			push(@{$tables_names->{$set}},$table);
			$documentation->{$set}{'tables'}{$table} = { 'desc' => '', 'column' => [], 'see' => [], 'info' => [] };
			$tag = $tag_content = '';		
		}
		# Description (used for both set, table and info tags)
		elsif ($doc =~ /^\@(desc)\s*(.+)$/i) {
			fill_documentation ($1,$2);
		}
		# Column
		elsif ($doc =~ /^\@(column)\s*(.+)$/i) {
			fill_documentation ($1,$2);
		}
		# See other tables
		elsif ($doc =~ /^\@(see)\s*(.+)$/i) {
			fill_documentation ($1,$2);	
		}
		# Addtional information block
		elsif ($doc =~ /^\@(info)\s*(.+)$/i) {
			fill_documentation ();
			$info = $2;
			next;
		}
		# End of documentation
		elsif ($doc =~ /^\*\//) { # End of the documentation block
			fill_documentation (); # Add the last tag content to the documentation hash
			$in_doc=0;
			next; 
		}
		elsif ($doc =~ /^\s*(.+)$/) { # If a tag content is split in several lines
			$tag_content .= " $1";
		}
	}
	
	## Parsing of the SQL table to fetch the columns types ##
	elsif ($in_table==1) {
	
		## INDEXES ##
		if ($doc =~ /^\s*(primary\skey)\s*\((.+)\)/i or $doc =~ /^\s*(unique)\s*\((.+)\)/i){ # Primary or unique
			add_column_index($1,$2);
			next;
		}
		elsif ($doc =~ /^\s*(unique\s)?(key)\s(\w+)\s*\((.+)\)/i) { # Keys and indexes
			add_column_index("$1$2",$4,$3);
			next;
		}
		
		## TYPES & DEFAULT VALUES ##
		my $col_name = '';
		my $col_type = '';
		my $col_def  = '';
		
		# All the type is contained in the same line (type followed by parenthesis)
		if ($doc =~ /^\s*(\w+)\s+(\w+\s?\(.*\))/ ){
			$col_name = $1;
			$col_type = $2;
			if ($doc =~ /default\s*([^,\s]+)\s*.*(,|#).*/i) { $col_def = $1; } # Default value
			add_column_type_and_default_value($col_name,$col_type,$col_def);
		}
		
		# The type is written in several lines
		elsif ($doc =~ /^\s*(\w+)\s+(enum|set)(\s?\(.*)/i){ # The content is split in several lines
			$col_name=$1;
			$col_type="$2$3<br />";
			my $end_type = 0;
			while ($end_type != 1){
				my $line = <SQLFILE>;
				chomp $line;
				
				if ($line =~ /\)/) { # Close parenthesis
					$end_type=1; 
					$line =~ /^\s*(.+)\)/;
					$col_type .= "$1)"; 
				}
				else { # Add the content of the line
					$line =~ /^\s*(.+)/;
					$col_type .= $1.'<br />';
				}
				if ($line =~ /default\s*([^,\s]+)\s*.*(,|#).*/i) { $col_def = $1; } # Default value
			}
			add_column_type_and_default_value($col_name,$col_type,$col_def);
		}
		
		# All the type is contained in the same line (type without parenthesis)
		elsif ($doc =~ /^\s*(\w+)\s+(\w+)/ ){
			$col_name = $1;
			$col_type = $2;
			if ($doc =~ /default\s*([^,\s]+)\s*.*(,|#).*/i) { $col_def = $1;} # Default value
			add_column_type_and_default_value($col_name,$col_type,$col_def);
		}
		
		elsif ($doc =~ /\).*;/) { # End of the sql table definition
			if (scalar @{$documentation->{$set}{'tables'}{$table}{column}} > $count_sql_col) {
				print STDERR "Description of a non existant column in the table $table!\n";
			}
			$in_table=0;
			$count_sql_col = 0;
			$table='';
		}
	}
}
close(SQLFILE);


############
### Core  ##
############

# Sort the tables names by alphabetic order
if ($sort_tables == 1) {
	while ( my($set_name,$tables) = each (%{$tables_names})) {
		@{$tables} = sort(@{$tables});
	}
}

# List of tables by set of tables
my $html_content = display_tables_list();
my $table_count = 1;
my $col_count = 1;

while ( my ($set_name,$tables) = each (%{$tables_names})) {
	# Set display	
	if ($set_flag == 1 and $set_name ne 'default') {
		$html_content .= qq{<hr />\n<h2>Set $set_name</h2>\n};
		my $set_desc = $documentation->{$set_name}{'desc'};		
		$html_content .= qq{<p>$set_desc</p>} if (defined($set_desc));
	}
	if ($set_name eq 'default' and defined($documentation->{$set_name}{'info'})) {
		$html_content .= qq{<h2>Additional information about the schema</h2>\n};
	}	
	$html_content .= add_info($documentation->{$set_name}{'info'});
	
	# Tables display
	foreach my $t_name (@{$tables}) {
		$html_content .= add_table_name($t_name);
		$html_content .= add_description($documentation->{$set_name}{'tables'}{$t_name}{desc});
		$html_content .= add_info($documentation->{$set_name}{'tables'}{$t_name}{info});	
		$html_content .= add_columns($t_name,@{$documentation->{$set_name}{'tables'}{$t_name}{column}});
		$html_content .= add_see(@{$documentation->{$set_name}{'tables'}{$t_name}{see}});
	}
}	


## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $header."\n";
print HTML $html_content."\n";
print HTML $footer."\n";
close(HTML);




###############
##  Methods  ##
###############

sub display_tables_list {
	my $html = qq{<p>List of the tables:</p>\n};
	while ( my ($set_name,$tables) = each (%{$tables_names})) {
		my $count = scalar @{$tables};
		my $nb_col = ceil($count/$nb_by_col);
		my $table_count = 0;
		my $col_count = 1;
	
		if ($nb_col>3) { 
			while ($nb_col>3) {
				$nb_by_col += 5;
				$nb_col = ceil($count/$nb_by_col);
			}
			$nb_col = 3;
		}

		if ($set_flag == 1 and $set_name ne 'default') {
			$html .= qq{<h2>Set $set_name</h2>\n};
		}
		$html .= qq{<table><tr><td>\n <ul>\n};

		foreach my $t_name (@{$tables}) {
			if ($table_count == $nb_by_col and $col_count<$nb_col and $nb_col>1){
				$html .= qq{	</ul></td><td><ul>\n};
				$table_count = 0;
			}
			$html .= add_table_name_to_list($t_name);
			$table_count ++;
		}
		$html .= qq{		</ul>\n</td></tr></table>\n};
	}
	$html .= qq{<hr />\n};
	return $html;
}


# If the line starts by a @<tag>, the previous tag content is added to the documentation hash.
# This method allows to describe the content of a tag in several lines.
sub fill_documentation {
	my $t1 = shift;
	my $t2 = shift;
	
	if ($tag ne '') {
		# Description tag (info, table or set)
		if ($tag eq 'desc') {
			# Additional description
			if ($info ne '') {
				$tag_content = $info.'@info@'.$tag_content;
				# Table: additional information				
				if ($table ne '') {
					push(@{$documentation->{$set}{'tables'}{$table}{'info'}},$tag_content);
				}
				# Set of tables: additional information
				else {
					if (!$documentation->{$set}{'info'}) {
						$documentation->{$set}{'info'} = [];
					}
					push(@{$documentation->{$set}{'info'}},$tag_content);
				}
				$info = '';
			}
			# Table description
			elsif ($table ne '') {
				$documentation->{$set}{'tables'}{$table}{$tag} = $tag_content;
			}
			# Set description
			else {
				$documentation->{$set}{'desc'} = $tag_content;
			}
		}
		elsif ($tag eq 'column') {
			$tag_content =~ /(\w+)[\s\t]+(.+)/;
			
			my $column = { 'name'    => $1,
								'type'    => '',
			               'default' => '',
								'index'   => '',
						      'desc'    => $2
							 };
			push(@{$documentation->{$set}{'tables'}{$table}{$tag}},$column);
		}
		else{
			push(@{$documentation->{$set}{'tables'}{$table}{$tag}},$tag_content);
		}
	}
	# New tag initialized
	if ($t1) {
		$tag = $t1;
		$tag_content = $t2;
	}
	else {
		$tag = $tag_content = '';	
	}
}
 

sub add_table_name_to_list {
	my $t_name = shift;
	
	my $html = qq{		<li><a href="#$t_name"><b>$t_name</b></a></li>\n};
	return $html;
}

sub add_table_name {
	my $t_name = shift;
	
	my $html = qq{\n<br />\n<table style="border: 2px groove #CCCCCC;height:10px;background-color:#FAFAFF"><tr style="vertical-align:middle;height:10px">
<td style="width:500px;text-align:left;height:10px"><span id="$t_name" style="font-size:11pt;font-weight:bold">$t_name</span></td>
<td style="width:100px;text-align:right"><a id="a_$t_name" style="cursor:pointer;text-decoration:underline" onclick="show_hide('$t_name')">Show</a> columns</td>
</tr></table>\n};
	
	return $html;
}


sub add_description {
	my $desc = shift;
	return qq{<p>$desc<\p>\n};
}

sub add_info {
	my $infos = shift;
	my $html = '';
	
	foreach my $inf (@{$infos}) {
		my ($title,$content) = split('@info@', $inf);
		$html .= qq{	<table>
		<tr><td class="bg3">$title</td><td class="bg1">$content</td></tr>
	</table>\n};
	}
	return $html;
}

sub add_columns {
	my @cols = @_;
	my $table = shift @cols;
	my $display_style = $display_col{$display};
	
	my $html = qq{\n	<div id="div_$table" style="display:$display_style">\n
	<table style="border:1px outset #222222">
		<tr class="bg3 center"><th style="width:180px">Column</th><th style="width:150px">Type</th><th style="width:100px">Default value</th><th style="width:400px">Description</th><th style="width:150px">Index</th></tr>\n};
	my $bg = 1;
	foreach my $col (@cols) {
		my $name    = $col->{name};
		my $type    = $col->{type};
		my $default = $col->{default};
		my $desc    = $col->{desc};
		my $index   = $col->{index};
		$desc =~ /\@link\s?(\w+)/;
		my $table_to_link = qq{<a href="#$1">$1</a>};
		$desc =~ s/\@link\s?\w+/$table_to_link/;
		
		#$col =~ /^\s*(\w+)[\s\t]+(.+)\t+(.+)\t(.*)/;
		$html .= qq{		<tr class="bg$bg"><td><b>$name</b></td><td>$type</td><td>$default</td><td>$desc</td><td>$index</td></tr>\n};
		if ($bg==1) { $bg=2; }
		else { $bg=1; }
	}
	$html .= qq {</table>\n</div>\n};
	
	return $html;
}

sub add_see {
	my @sees = @_;
	my $html = '';
	
	if (scalar @sees) {
		$html .= qq{<p><b>See also:</b></p>\n<ul>\n};
		foreach my $see (@sees) {
			$html .= qq{	<li><a href="#$see">$see</a></li>\n};
		}
		$html .= qq{</ul>\n};
	}
	
	return $html;
}


sub add_column_index {
	my $idx_type = shift;
	my $idx_col  = shift;
	my $idx_name = shift;
	
	my $index = $idx_type;
	if (defined($idx_name)) {
		$index .= ": $idx_name";
	}
	
	my @idx_cols = split(',',$idx_col); # The index can involve several columns
	
	my %is_found = ();
	foreach my $i_col (@idx_cols) {
		$i_col =~ s/^\s+//; # Remove white spaces
		$i_col =~ s/\s+$//;
		
		$is_found{$i_col} = 0;
		foreach my $col (@{$documentation->{$set}{tables}{$table}{column}}) {
			if ($col->{name} eq $i_col) {
				if ($col->{index} ne '') {
					$col->{index} .= '<br />';
				}
				$col->{index} .= lc($index);
				$is_found{$i_col} = 1;
				last;
			}
		}
	}
	# Description missing
	while (my ($k,$v) = each(%is_found)) {
		if ($v==0) {
			print STDERR "The description of the column '$k' is missing in the table $table!\n";
		}
	}
}


sub add_column_type_and_default_value {
	my $c_name    = shift;
	my $c_type    = shift;
	my $c_default = shift;
	$count_sql_col ++;
	
	my $is_found = 0;
	foreach my $col (@{$documentation->{$set}{'tables'}{$table}{column}}) {
		if ($col->{name} eq $c_name) {
			#$col .= "\t$c_type\t";
			#$col .= "$c_default" if ($c_default ne ''); # Add the default value
			$col->{type} = $c_type;
			$col->{default} = $c_default if ($c_default ne '');
			$is_found = 1;
			last;
		}
	}
	# Description missing
	if ($is_found==0) {
		print STDERR "The description of the column '$c_name' is missing in the table $table!\n";
	}
}


sub usage {
	
  print qq{
  Usage: perl sql2html.pl [OPTION]
  
  Convert the SQL documentation into an HTML document.
	
  Options:

    -help		Print this message
      
    An input file must be specified. This file must be a SQL file, with the "Java-doc like documentation". 
    For more information, please visit the following page: 
    http://www.ebi.ac.uk/seqdb/confluence/display/EV/SQL+documentation

    -i          A SQL file name (Required)
    -o          An HTML output file name (Required) 
    -show_set   A flag to take into account the set of tables (1) or not (0). 
                By default, the value is set to 1.
    -sort       A flag to sort (1) or not (0) the tables by alphabetic order.
                By default, the value is set to 1.
  } . "\n";
  exit(0);
}
