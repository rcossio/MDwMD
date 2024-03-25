const express = require('express');
const router = express.Router();
const Bottleneck = require('bottleneck');
const xml2js = require('xml2js').parseStringPromise; 
const Diffusion = require('../models/diffusion'); 
const Reference = require('../models/reference'); 
const Protein = require('../models/protein');

// Initialize a new bottleneck limiter
const limiter = new Bottleneck({
  maxConcurrent: 1, 
  minTime: 500
});

// Define your routes here
router.get('/', (req, res) => {
  res.send('Hello, World!');
});

// Serve HTML Page
router.get('/', (req, res) => {
  res.sendFile(path.join(__dirname, '/public/index.html'));
});

// Handle POST Request
router.post('/manage-data', async (req, res) => {
  const { data, replacedEntry } = req.body;

  try {
    const newEntryId = await Diffusion.create({ ...data, active: true, timestamp: new Date().toISOString() });
    if (replacedEntry) {
      await Diffusion.updateOne({ _id: replacedEntry }, { $set: { supersededBy: newEntryId, active: false } });
    }
    res.send({ status: 'success', payload: 'Data saved to MongoDB' });
  } catch (err) {
    console.error('Error saving data:', err);
    res.status(500).send('Error saving data');
  }
});

// Handle PUT Request
router.put('/manage-data', async (req, res) => {
  const { id, data } = req.body;
  Diffusion.updateOne({ _id: id }, { $set: data })
    .then(() => res.send({ status: 'success', payload: 'Data updated in MongoDB' }))
    .catch(err => {
      console.error('Error updating data:', err);
      res.status(500).send('Error updating data');
    });
});

// Handle POST Request for Search
router.post('/search', async (req, res) => {
  const { query, species, referenceType } = req.body;

  let queryConditions = [{ active: true }, {
    $or: [
      { proteinName: new RegExp(query, 'i') },
      { accessionNumber: new RegExp(query, 'i') },
      { referenceId: new RegExp(query, 'i') }
    ]
  }];

  if (species) {
    const speciesList = species.split(',').map(item => item.trim()).filter(item => item !== '');
    if (speciesList.length > 0) {
      queryConditions.push({ species: { $in: speciesList } });
    }
  }

  if (referenceType) {
    queryConditions.push({ referenceType });
  }

  try {
    const data = await Diffusion.find({ $and: queryConditions }).sort({ accessionNumber: 1 });

    const publicDataTasks = data.map(async (item) => {
      let year;
      const reference = await Reference.findOne({ $or: [{ pubmedId: item.referenceId }, { DOI: item.referenceId }] });

      if (reference) {
          year = reference.year;
      } else {
          year = await fetchPublicationYear(item.referenceId, item.referenceIdType);

          // TO FIX: Save the reference data to the database WHEN POSTING THE DATA, NOT WHEN SEARCHING
          try {
            if (item.referenceIdType === 'pmid') {
              await Reference.create({ pubmedId: item.referenceId, year: year });
            } else if (item.referenceIdType === 'doi') {
              await Reference.create({ DOI: item.referenceId, year: year });
            }
          } catch (error) {
            console.error('Error saving reference data:', error);
          }
      }

      return {
        _id: item._id,
        diffusionCoefficient: item.diffusionCoefficient,
        diffusionError: item.diffusionError,
        diffusionUnit: item.diffusionUnit,
        manualProcessing: item.manualProcessing,
        otherManualProcessing: item.otherManualProcessing,
        method: item.method,
        otherMethod: item.otherMethod,
        proteinName: item.proteinName,
        species: item.species,
        accessionNumber: item.accessionNumber,
        referenceId: item.referenceId,
        referenceIdType: item.referenceIdType,
        referenceType: item.referenceType,
        comment: item.comment,
        pdbStructures: item.pdbStructures,
        discardedPdbStructures: item.discardedPdbStructures,
        year: year
      };
    });

    const publicData = await Promise.all(publicDataTasks);

    const uniqueSpecies = [...new Set(publicData.map(item => item.species))].map(Number).sort((a, b) => a - b);
    const uniqueReferenceTypes = [...new Set(publicData.map(item => item.referenceType))].sort((a, b) => a.localeCompare(b));

    res.json({
      status: 'success',
      payload: {
        results: publicData, 
        species: uniqueSpecies,
        referenceTypes: uniqueReferenceTypes
      }
    });
  } catch (error) {
    console.error('Error during search:', error);
    res.status(500).send('Error during search');
  }
});


async function fetchPublicationYear(refId, refIdType) {
  return limiter.schedule(async () => {
    let url = refIdType === 'pmid' ? 
              `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=${refId}&retmode=xml` : 
              `https://api.crossref.org/works/${refId}`;

    try {
      const response = await fetch(url);
      const text = await response.text();

      if (refIdType === 'pmid') {
        const result = await xml2js(text);
        const year = result.PubmedArticleSet.PubmedArticle[0].MedlineCitation[0].Article[0].Journal[0].JournalIssue[0].PubDate[0].Year[0];
        return year;
      } else if (refIdType === 'doi') {
        const json = JSON.parse(text);
        const year = json.message['published-print']?.['date-parts'][0][0] || json.message['published-online']?.['date-parts'][0][0];
        return year;
      }
    } catch (error) {
      console.error('Fetch error:', error);
      throw new Error('Failed to fetch publication year');
    }
  });
}

// Handle GET Request to get by id
router.get('/find-by-id/:id', (req, res) => {
  const { id } = req.params;

  Diffusion.findById(id)
    .then(data => res.json({
      status: 'success',
      payload: data
    }))
    .catch(err => {
      console.error('Error during findById:', err);
      res.status(500).send('Error during findById');
    });
});
  

// Handle POST Request
router.post('/protein', async (req, res) => {
  const { data } = req.body;

  try {
    await Protein.create({ ...data, active: true, timestamp: new Date().toISOString() });

    res.send({ status: 'success', payload: 'Data saved to MongoDB' });
  } catch (err) {
    console.error('Error saving data:', err);
    res.status(500).send('Error saving data');
  }
});

module.exports = router;
