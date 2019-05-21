# v0.1

## 5-17-19
- Add the position of the C>U change in the tsv file
- Remove CUAC, CCUC, etc. requirement for 'best' loop. Now just require a C in the 0 position and a C or U in the -1 position.

# v0.2

## 5-21-19
- Add option to use best tetraloop from the 2015 Bora and Shraddha paper https://www.nature.com/articles/ncomms7881.
- Add function to check for bulge at the -2 position automatically. If true, it will continue the stem search.
- Change energy and stem ranknig to ascending and descending, respectively.
- Add larger subsequence energy estimates
	- Commented for now. 
	- Review and test later

