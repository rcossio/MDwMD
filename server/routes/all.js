const express = require('express');
const router = express.Router();
const Bottleneck = require('bottleneck');
const xml2js = require('xml2js');
const DiffusionData = require('../dao/models'); //TO DO: migrate to a data access layer

// Initialize a new bottleneck limiter
const limiter = new Bottleneck({
  maxConcurrent: 1, // Maximum number of jobs running at the same time
  minTime: 500 // Minimum time (in milliseconds) between job executions
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
    const newEntryId = await DiffusionData.create({ ...data, active: true, timestamp: new Date().toISOString() });
    if (replacedEntry) {
      await DiffusionData.updateOne({ _id: replacedEntry }, { $set: { supersededBy: newEntryId, active: false } });
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
  DiffusionData.updateOne({ _id: id }, { $set: data })
    .then(() => res.send({ status: 'success', payload: 'Data updated in MongoDB' }))
    .catch(err => {
      console.error('Error updating data:', err);
      res.status(500).send('Error updating data');
    });
});

// Handle GET Request for Search
router.post('/search', (req, res) => {
  const { query, species, referenceType } = req.body;

  let queryConditions = [
    { active: true },
    {
      $or: [
        { proteinName: new RegExp(query, 'i') },
        { accessionNumber: new RegExp(query, 'i') },
        { referenceId: new RegExp(query, 'i') }
      ]
    }
  ];

  // Add species filter only if species is not empty and not null
  const speciesList = species.split(',').map(item => item.trim()).filter(item => item !== '');

  if (speciesList && speciesList.length > 0) {
    queryConditions.push({ species: { $in: speciesList } });
  }

  // Add referenceType filter only if referenceType is not null
  if (referenceType) {
    queryConditions.push({ referenceType });
  }

  DiffusionData.find({ $and: queryConditions })
    .sort({ accessionNumber: 1 })
    .then(data => {
      // Extract species from the results and create a unique list
      const species = new Set(data.map(item => item.species));
      const speciesList = Array.from(species).map(Number).sort((a, b) => a - b);

      // Extract referenceType from the results and create a unique list
      const referenceTypes = new Set(data.map(item => item.referenceType));
      const referenceTypesList = Array.from(referenceTypes).sort((a, b) => a.localeCompare(b));

      // Send the response with both results and species list
      res.json({
        status: 'success',
        payload: {
          results: data,
          species: speciesList,
          referenceTypes: referenceTypesList
        }
      });
    })
    .catch(err => {
      console.error('Error during search:', err);
      res.status(500).send('Error during search');
    });
});

// Handle GET Request to get by id
router.get('/find-by-id/:id', (req, res) => {
  const { id } = req.params;

  DiffusionData.findById(id)
    .then(data => res.json({
      status: 'success',
      payload: data
    }))
    .catch(err => {
      console.error('Error during findById:', err);
      res.status(500).send('Error during findById');
    });
});

// Handle GET Request for fetching publication year
router.get('/publication-year', (req, res) => {
    const { refId, refType } = req.query;
  
    const fetchPublicationYear = async () => {
      let url, response, text;
      if (refType === 'pmid') {
        url = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=${refId}&retmode=xml`;
      } else if (refType === 'doi') {
        url = `https://api.crossref.org/works/${refId}`;
      }
      
      try {
        const fetch = (await import('node-fetch')).default;
        response = await fetch(url);
        text = await response.text();
      } catch (error) {
        console.error('Fetch error:', error);
        return Promise.reject('Failed to fetch data');
      }
  
      if (refType === 'pmid') {
        return new Promise((resolve, reject) => {
          xml2js.parseString(text, (err, result) => {
            if (err) {
              reject('Error parsing XML');
            } else {
              try {
                const year = result.PubmedArticleSet.PubmedArticle[0].MedlineCitation[0].Article[0].Journal[0].JournalIssue[0].PubDate[0].Year[0];
                resolve(year);
              } catch (parseError) {
                reject('Error extracting publication year');
              }
            }
          });
        });
      } else if (refType === 'doi') {
        try {
          const json = JSON.parse(text);
          const year = json.message['published-print']?.['date-parts'][0][0] || json.message['published-online']?.['date-parts'][0][0];
          return Promise.resolve(year);
        } catch (parseError) {
          console.error('Error parsing JSON:', parseError);
          return Promise.reject('Error extracting publication year from DOI');
        }
      }
    };
  
    limiter.schedule(fetchPublicationYear)
      .then(year => res.json({ year }))
      .catch(error => {
        console.error(error);
        res.status(500).send(error);
    });
});
  
  

module.exports = router;
