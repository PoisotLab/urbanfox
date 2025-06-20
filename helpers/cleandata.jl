using SpeciesDistributionToolkit
import CSV

# Read the line data from Pharos export
raw = CSV.File("data/pharos_data.csv")

# Now we loop and create an occurrence from each record
raw[2]

function depharosizer(record)
    return Occurrence(
        what = record.detection_target,
        where = (record.longitude, record.latitude),
        when = record.collection_date,
        presence = record.detection_outcome == "positive"
    )
end

outcomes = Occurrences(depharosizer.(raw))
OccurrencesInterface.save("data/occurrences.json", outcomes)