#!/bin/sh

# Run after each job in this project has run.
# Tracks usage
# Checks for the usage file in /home/ekoch/sims/

# Must be executable:
# chmod u+x epilogue.script
# chmod go-rwx epilogue.script

echo "Running epilogue."
echo "Epilogue Args:"
echo "Job ID: $1"
echo "User ID: $2"
echo "Group ID: $3"
echo "Job Name: $4"
echo "Session ID: $5"
echo "Resource List: $6"
echo "Resources Used: $7"
echo "Queue Name: $8"
echo "Account String: $9"
echo ""

echo "Appending results to usage file in /home/ekoch/sims/astrostat_usage.log"

if [ ! -f /home/ekoch/astrostat_usage.log ]; then
    # Need to create the file.
    echo "Job ID, User ID, Group ID, Job Name, Session ID, Resource List, Resource Used, Queue Name, Account String" >> /home/ekoch/sims/astrostat_usage.log
fi

echo "$1, $2, $3, $4, $5, $6, $7, $8, $9" >> /home/ekoch/sims/astrostat_usage.log

exit 0