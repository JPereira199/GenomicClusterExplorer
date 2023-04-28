#!/usr/bin/perl

# Especie de fastacomand. si le das file_out, te imprime un multifasta con todos los de la lista.
# Si no le das file_out, te imprime cada uno de los que encuentra en un fasta individual con su nombre.
# En ambos casos imprime (en pantalla) los nombres que no encontró.
# En el fasta de entrada, los nombres de las secuencias deben estar separados por uno o más espacios de la descripción
# (siempre hablando de la línea que empieza con >), si es que hay descripción.

use strict;

my $usage = "\nperl $0 fasta_in lista_nombres [fasta_out]\n"; 
if (scalar @ARGV == 0){print $usage};
my $flag = 0;

my $fasta_in = $ARGV[0];
my %hash_nc_sec = leer_multi_fasta($fasta_in);
my %hash_nc_nl;


# print "primer hash\n";
# while (my($key,$value) = each(%hash_nc_sec)){
	# print "$key........$value\n";
# }
# print "segundo hash\n";
# while (my($key,$value) = each(%hash_nc_nl)){
	# print "$key--------$value\n";
# }



if ($ARGV[2]){
	print "No encontrados: ";
	open (OUT,">$ARGV[2]");
	open (IN, $ARGV[1]) or die "No existe el archivo $ARGV[1], o falta especificar un segundo argumento\n";
	while(<IN>){
		chomp $_;
		if(-exists $hash_nc_sec{$_}){
			print OUT ">$hash_nc_nl{$_}\n$hash_nc_sec{$_}\n";	
		}
		else {
			print "\n$_";
			$flag = 1;
		}
	}
	close OUT;
}
else{
	print "No encontrados: ";
	open (IN, $ARGV[1]) or die "No existe el archivo $ARGV[1], o falta especificar un segundo argumento\n";
	while(<IN>){
		chomp $_;
		if(-exists $hash_nc_sec{$_}){
			open (OUT,">$_");
			print OUT ">$hash_nc_nl{$_}\n$hash_nc_sec{$_}\n";
			close OUT;
		}
		else {
			print "\n$_";
			$flag = 1;
		}
	}
}
if($flag == 0){ print "0\n"}
else {print "\n"};

####################################################################################################################

sub leer_multi_fasta {
	my $archivo_multi_fasta = $_[0];
	my %hash = ();
	my $secuencia = '';
	my $name_sec = '';
	my $name_corto = '';
	
	open (READ,$archivo_multi_fasta) or die "No existe el archivo $archivo_multi_fasta\n";
	while(<READ>){
		chomp ($_);
		
		if(/>(.*)/) {
        	if($secuencia) {
        		$hash{$name_corto} = $secuencia;
        		$hash_nc_nl{$name_corto} = $name_sec;
        	}
        	$name_sec = $1;		#	print $name_sec,"\n"; 
        	$name_sec =~ /^(\S+)/;
        	$name_corto = $1;	#	pint "$name_corto\n";
			$secuencia = '';
		}
		else {
			$_ = uc $_;
			$_ =~ s/[0-9]//g;
			$_ =~ s/\s//g;
        	$secuencia .= $_;
   		}
	}
	$hash{$name_corto} = $secuencia;
	$hash_nc_nl{$name_corto} = $name_sec;
	close(READ);
	%hash;
}
