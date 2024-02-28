const express = require('express');
const bodyParser = require('body-parser');
const mongoose = require('mongoose');
const xml2js = require('xml2js');
const Bottleneck = require('bottleneck');
require('dotenv').config();

const app = express();
const PORT = 3000;

// Initialize a new bottleneck limiter
const limiter = new Bottleneck({
  maxConcurrent: 1, // Maximum number of jobs running at the same time
  minTime: 500 // Minimum time (in milliseconds) between job executions
});

// MongoDB Schema and Model setup
const diffusionSchema = new mongoose.Schema({
  diffusionCoefficient: Number,
  diffusionError: Number,
  diffusionUnit: String,
  manualProcessing: String,
  otherManualProcessing: String,  
  method: String,
  otherMethod: String,
  proteinName: String,
  species: String,
  accessionNumber: String,
  pubmedId: String,
  referenceType: String,
  comment: String, 
  timestamp: String,
  userName: String 
});

// Create Mongoose Model
const DiffusionData = mongoose.model('data_from_experiments', diffusionSchema);

app.use(bodyParser.json());
app.use(express.static('public'));

// Serve HTML Page
app.get('/', (req, res) => {
  res.sendFile(path.join(__dirname, '/public/index.html'));
});

// Handle POST Request
app.post('/update-data', (req, res) => {
  const newData = req.body;

  DiffusionData.create(newData) 
  .then(() => {
    res.send({ status: 'success', payload: 'Data saved to MongoDB' });
  })
  .catch((err) => {
    console.error('Error saving data:', err);
    res.status(500).send('Error saving data');
  });
});

// Handle GET Request for Search
app.get('/search', (req, res) => {
  const { query } = req.query;

  DiffusionData.find({
    $or: [
      { proteinName: new RegExp(query, 'i') },
      { accessionNumber: new RegExp(query, 'i') }
    ]
  })
  .sort({ accessionNumber: 1 })
  .then(data => {
    res.send({ status: 'success', payload: data });
  })
  .catch(err => {
    console.error('Error during search:', err);
    res.status(500).send('Error during search');
  });
});

// Handle GET Request for fetching publication year
app.get('/publication-year', (req, res) => {
  const { pmid } = req.query;

  const fetchPublicationYear = async () => {
    const url = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=${pmid}&idtype=pmid&retmode=xml`;

    // Fetch response using node-fetch
    const response = await (await (await import('node-fetch')).default(url)).text();
    return new Promise((resolve, reject) => {
      xml2js.parseString(response, (err, result) => {
        if (err) {
          reject('Error parsing XML');
          console.error(err);
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
  };

  limiter.schedule(fetchPublicationYear)
    .then(year => res.json({ year }))
    .catch(error => {
      console.error(error);
      res.status(500).send(error);
    });
});

mongoose.connect(process.env.MONGODB_URI)
.then(() => console.log('Connected to MongoDB'))
.catch(err => console.error('Error connecting to MongoDB:', err));

app.listen(PORT, () => {
  console.log(`Server running on http://localhost:${PORT}`);
});
