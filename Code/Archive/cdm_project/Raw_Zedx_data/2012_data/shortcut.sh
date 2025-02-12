#!/bin/bash
echo "This is a folder file manipulation script" 

ls -1F > mylisting.txt
sed 's/\'_wxdata.csv'//' mylisting.txt > new.txt
sed 's/_/ /' new.txt > new_listing.txt
sed 's/_/ /' new_listing.txt > new_listing2.txt
sed 's/ /,/' new_listing2.txt > new_listing.csv
