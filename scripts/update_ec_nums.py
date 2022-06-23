from brendapy import BrendaParser


# Reload the BRENDA Parser to get the latest EC numbers
# BRENDA_PARSER = BrendaParser()


# keys = random.sample(BRENDA_PARSER.keys(), 2)

keys = ['2.4.3.1', '4.5.6.7']

print ('called')

print (snakemake.output)

with open (snakemake.output) as ec_nums_txt:
	ec_nums_txt.write (f'ec_nums = [{keys}]')
